import fwdpy11

import fwdpy11.genetic_values
import fwdpy11.model_params
import fwdpy11.wright_fisher
import fwdpy11.multilocus
import fwdpy11.sampling
import numpy as np
import pandas as pd
import sqlite3
import os
import sys
import argparse
import scipy.stats
import scipy.optimize
from collections import namedtuple
import concurrent.futures
import libsequence.variant_matrix as vmatrix
import libsequence.summstats as sstats

# The next two get imported for testing only
# import libsequence.polytable
# import libsequence.windows

# A streamlined 10-locus simulator.
# "plarge" is an input parameter.
# We do a one-stop-shop in terms of what we
# track.  No pickling or whatnot.
#
# Currently requires dev branches of all tools.


def make_parser():
    parser = argparse.ArgumentParser(
        description="Simple multi-locus simulation of single optimum shift")

    parser.add_argument('--N', '-N', type=int, default=5000,
                        help="Diploid population size")
    parser.add_argument('--mu', '-m', type=float, default=None,
                        help="mutation rate to mutations "
                        "affecting trait value")
    parser.add_argument('--theta', '-t', type=float,
                        default=100.0, metavar="THETA", help="4Nu")
    parser.add_argument('--rho', '-r', type=float,
                        default=100.0, metavar="RHO", help="4Nr")
    parser.add_argument('--seed', '-S', type=int,
                        default=None, help="RNG seed")
    parser.add_argument('--nsam', type=int, default=None,
                        help="Number of diploids to sample for genome scan")
    parser.add_argument('--ncores', '-n', type=int, default=1)
    parser.add_argument('--nreps', '-nr', type=int, default=1)
    parser.add_argument('--nloci', '-nl', type=int, default=10)
    parser.add_argument('--plarge', '-si', type=float, default=None)
    parser.add_argument('--opt', '-o', type=float, default=0.0)
    parser.add_argument('--qoutfile', type=str,
                        help="Output file name for quantitative trait data")
    parser.add_argument('--goutfile', type=str,
                        help="Output file name for genome scan data")
    parser.add_argument('--foutfile', type=str,
                        help="Output file name for recent fixation data")
    parser.add_argument('--gamma', type=float, default=None,
                        help="Use Gamma DES with specific shape.")
    return parser


def gamma_hat(VS, mu):
    return 2.0 * np.sqrt(2.0) * np.sqrt(VS * mu)


def generate_gaussian_function_to_minimize(ghat, plarge):
    return lambda x: abs(2.0 *
                         (1.0 - scipy.stats.norm.cdf(ghat, scale=x)) - plarge)


def get_gaussian_sigma(F):
    res = scipy.optimize.minimize_scalar(F, bounds=(0, 100), method='bounded')
    return res.x


def minimize_gamma_cdf(x, shape, ghat, plarge):
    return (1.0 - scipy.stats.gamma.cdf(ghat, a=shape, scale=x / shape) - plarge)**2


# namedtuple('QData', ['generation', 'vg', 'zbar', 'g_per_locus'])
QData = None
GSData = namedtuple(
    'GSData', ['generation', 'locus', 'window', 'thetapi',
               'tajd', 'hprime', 'H1', 'H12', 'H2H1', 'meanz'])
SGV = namedtuple('SGV', ['s', 'g', 'pos', 'nhaps_at_shift'])


class Sampler(object):
    def __init__(self, rng, nsam):
        self.qdata = []
        self.gsdata = []
        self.sgv = []
        self.rng = rng
        self.nsam = nsam

    def _gvalue_per_locus(self, pop):
        """
        Strict additive!!!
        """
        rv = np.zeros(len(pop.locus_boundaries))
        mutations = np.array(pop.mutations.array())
        for dip in pop.diploids:
            for i, j in zip(dip, range(len(pop.locus_boundaries))):
                for gi in [i.first, i.second]:
                    g = np.array(pop.gametes[gi].smutations, copy=False)
                    rv[j] += mutations['s'][g].sum()
        return rv / pop.N

    def _selected_mut_keys(self, pop):
        md = {i: 0 for i in range(len(pop.mcounts)) if
              pop.mcounts[i] > 0 and pop.mcounts[i] < 2 * pop.N and
              pop.mutations[i].neutral is False}
        if len(md) > 0:
            for g in pop.gametes:
                if g.n > 0:
                    for k in md.keys():
                        if k in g.smutations:
                            md[k] += 1
            for k, v in md.items():
                self.sgv.append(
                    SGV(pop.mutations[k].s,
                        pop.mutations[k].g,
                        pop.mutations[k].pos, v))

    def _merge(self, sample):
        npos = np.array(sample.neutral.positions)
        spos = np.array(sample.selected.positions)

        if len(spos) == 0 and len(npos) == 0:
            return None
        if len(spos) == 0 and len(npos) != 0:
            return vmatrix.VariantMatrix(np.array(sample.neutral,
                                                  copy=False), npos)
        if len(spos) != 0 and len(npos) == 0:
            return vmatrix.VariantMatrix(np.array(sample.selected,
                                                  copy=False), spos)

        # We need to merge the arrays into one glorious numpy mashup
        mpos = np.concatenate((npos, spos))
        data = np.concatenate((np.array(sample.neutral, copy=False),
                               np.array(sample.selected, copy=False)))
        mpos_idx = np.array(
            sorted([i for i in range(len(mpos))], key=lambda x: mpos[x]))
        return vmatrix.VariantMatrix(data[mpos_idx, :], mpos[mpos_idx])

    def _mean_zscore_by_window(self, sample, locus_boundaries):
        """
        Gets the mean z-score associated with nSL for all
        (locus, window) pairs.  The first/last window
        of each locus is used for binning.
        """
        # Get nSL values for each variant at each locus
        nSL_tips = None
        for lb in locus_boundaries:
            # Get the stats for the first and last window
            w = sample.window(*lb)
            w0 = sample.window(lb[0], lb[0] + 1)
            w1 = sample.window(lb[1] - 1, lb[1])
            nsl0 = sstats.nsl(w0, 0)
            nsl1 = sstats.nsl(w1, 0)
            if nSL_tips is None:
                nSL_tips = np.array(nsl0)
                nSL_tips = np.concatenate(
                    (nSL_tips, np.array(nsl1, copy=False)))
            else:
                nSL_tips = np.concatenate((nSL_tips,
                                           np.array(nsl0, copy=False),
                                           np.array(nsl1, copy=False)))

        valid_tip_stats = np.where(np.logical_and(
            np.isfinite(nSL_tips['nsl']), nSL_tips['core_count'] > 3))[0]
        tip_DAF = nSL_tips['core_count'][valid_tip_stats]
        tip_STAT = nSL_tips['nsl'][valid_tip_stats]
        bins = np.digitize(tip_DAF, np.arange(0, sample.nsam, 10))
        mean_sd = {}
        for i in np.unique(bins):
            w = np.where(bins == i)[0]
            if len(w) > 0:
                m = tip_STAT[w].mean()
                sdev = tip_STAT[w].std()
                if np.isfinite(sdev) == 1 and sdev > 0.0:
                    mean_sd[i] = (m, sdev)

        # Now, we need to return a mapping
        # of (locus, window) -> mean z score
        mapping = dict()
        for lb, locus in zip(locus_boundaries, range(len(locus_boundaries))):
            wdw = 1
            for l in range(int(lb[0]) + 1, int(lb[1]) - 1):
                window = sample.window(l, l + 1)
                nsl = np.array(sstats.nsl(window, 0))
                valid_tip_stats = np.where(np.logical_and(
                    np.isfinite(nsl['nsl']), nsl['core_count'] > 3))[0]
                stats = nsl['nsl'][valid_tip_stats]
                daf = nsl['core_count'][valid_tip_stats]
                bins = np.digitize(daf, np.arange(0, sample.nsam, 10))
                zscores_win = np.array([])
                for i in np.unique(bins):
                    if i in mean_sd:
                        w = np.where(bins == i)[0]
                        m = mean_sd[i][0]
                        sdev = mean_sd[i][1]
                        z = (stats[w] - m) / sdev
                        zscores_win = np.concatenate((zscores_win, z))
                if len(zscores_win) > 0:
                    mz = zscores_win.mean()
                else:
                    mz = np.nan
                mapping[(locus, wdw)] = mz
                wdw += 1
            # Window 0 and 10 get nan
            assert (locus, 0) not in mapping, "logic fail"
            assert (locus, 10) not in mapping, "logic fail"
            mapping[(locus, 0)] = np.nan
            mapping[(locus, 10)] = np.nan

        return mapping

    def __call__(self, pop):
        if pop.generation >= 9 * pop.N:
            c1 = pop.generation % 100 == 0
            c2 = pop.generation < 11 * pop.N and pop.generation % 10 == 0.0
            c3 = pop.generation > 10 * pop.N and pop.generation < 10 * pop.N + 100
            if pop.generation == 10 * pop.N:
                # This shift just happened,
                # so record keys for all non-neutral variants
                self._selected_mut_keys(pop)
            if c1 or c2 or c3:
                if pop.generation % 10 == 0:
                    # Do this less often b/c it is relatively slow
                    # to get the genetic values per locus
                    md = np.array(pop.diploid_metadata, copy=False)
                    gpl = self._gvalue_per_locus(pop)
                    self.qdata.append(
                        QData(pop.generation, md['g'].var(), md['g'].mean(), *(i for i in gpl)))
                sample = pop.sample(self.rng, self.nsam)
                # For testing, to make sure I don't eff up the windowing below :)
                # sample4libseq = fwdpy11.sampling.matrix_to_sample(sample)
                # nl = fwdpy11.sampling.separate_samples_by_loci(
                #     pop.locus_boundaries, sample4libseq[0])
                # sl = fwdpy11.sampling.separate_samples_by_loci(
                #     pop.locus_boundaries, sample4libseq[1])
                # test_results = {}
                # for neut, sel, lb, locus in zip(nl, sl, pop.locus_boundaries, range(len(nl))):
                #     neut.extend([i for i in sel])
                #     neut = sorted(neut, key=lambda x: x[0])
                #     sd = libsequence.polytable.SimData(neut)
                #     w = libsequence.windows.Windows(sd, 1, 1, *lb)
                #     I = 0
                #     for wi in w:
                #         ps = libsequence.summstats.PolySIM(wi)
                #         d = ps.tajimasd()
                #         hp = ps.hprime()
                #         test_results[(locus, I)] = (d, hp)
                #         I += 1

                sample = self._merge(sample)
                # assert np.all(
                #     sample.positions[:-1] <= sample.positions[1:]), "sort error"
                mean_z = self._mean_zscore_by_window(
                    sample, pop.locus_boundaries)
                # Go over loci:
                for lb, locus in zip(pop.locus_boundaries,
                                     range(len(pop.locus_boundaries))):
                    w = sample.window(*lb)
                    p = np.array(w.positions)
                    ac = w.count_alleles()
                    windex = 0
                    for left in np.arange(lb[0], lb[1], 1):
                        idx = np.where((p >= left) & (p < left + 1.0))[0]
                        D = sstats.tajd(ac[idx])
                        Hp = sstats.hprime(ac[idx], 0)
                        pi = sstats.thetapi(ac[idx])
                        subw = w.window(left, left + 1)
                        g = sstats.garud_statistics(subw)
                        self.gsdata.append(GSData(pop.generation,
                                                  locus, windex, pi, D, Hp,
                                                  g.H1, g.H12, g.H2H1,
                                                  mean_z[(locus, windex)]))
                        # print(locus, w, D, Hp, test_results[(locus, w)])
                        windex += 1


def runsim(args):
    args, repid, repseed = args
    rng = fwdpy11.GSLrng(repseed)
    locus_boundaries = [(float(i + i * 11), float(i + i * 11 + 11))
                        for i in range(args.nloci)]
    NANC = args.N
    # In revision, change from using gamma_hat
    # to (2Ngamma^2)/(2VS) >= 100.
    # Keep the variable name to minimize 
    # the refactoring...
    # ghat = gamma_hat(1.0, args.mu)
    ghat = np.sqrt(100./5000.)
    sregions = None
    if args.gamma is None:
        F = generate_gaussian_function_to_minimize(ghat, args.plarge)
        sigma_gamma = get_gaussian_sigma(F)
        sregions = [[fwdpy11.GaussianS(j[0] + 5., j[0] + 6.,
                                       args.mu, sigma_gamma, coupled=False)]
                    for i, j in zip(range(args.nloci), locus_boundaries)]
    else:
        # Get a shape param for the gamma
        res = scipy.optimize.minimize_scalar(minimize_gamma_cdf, bounds=(
            0, 100), method='bounded', args=(args.gamma, ghat, args.plarge))
        sregions = [[fwdpy11.GammaS(j[0] + 5., j[0] + 6.,
                                    args.mu, -1.0 * res.x, args.gamma, coupled=False),
                     fwdpy11.GammaS(j[0] + 5., j[0] + 6.,
                                    args.mu, res.x, args.gamma, coupled=False)]
                    for i, j in zip(range(args.nloci), locus_boundaries)]
    print(sregions)
    nregions = [[fwdpy11.Region(j[0], j[1],
                                args.theta / (4. * float(NANC)), coupled=True)]
                for i, j in zip(range(args.nloci), locus_boundaries)]
    recregions = [[fwdpy11.Region(j[0], j[1],
                                  args.rho / (4. * float(NANC)), coupled=True)]
                  for i, j in zip(range(args.nloci), locus_boundaries)]
    interlocus_rec = fwdpy11.multilocus.binomial_rec([0.5] * (args.nloci - 1))
    nlist = np.array([NANC] * 15 * NANC, dtype=np.uint32)
    env = [(0, 0, 1), (10 * NANC, args.opt, 1)]
    mutrates_s = [args.mu / float(args.nloci)] * args.nloci
    mutrates_n = [10 * args.theta / float(4 * NANC)] * args.nloci
    recrates = [10 * args.rho / float(4 * NANC)] * args.nloci
    gssmo = fwdpy11.genetic_values.GSSmo(env)
    gv = fwdpy11.genetic_values.MlocusAdditive(2.0, gssmo)
    pdict = {'nregions': nregions,
             'sregions': sregions,
             'recregions': recregions,
             'demography': nlist,
             'interlocus_rec': interlocus_rec,
             'gvalue': gv,
             'rates': (mutrates_n, mutrates_s, recrates),
             'prune_selected': False
             }
    params = fwdpy11.model_params.ModelParams(**pdict)
    pop = fwdpy11.MlocusPop(args.nloci, locus_boundaries)
    sampler = Sampler(rng, args.nsam)
    fwdpy11.wright_fisher.evolve(rng, pop, params, sampler)
    sampler.rng = None  # Cannot be pickled and is not needed
    return repid, pop, sampler


def classify_fixations(pop, sgv, mu):
    sgvdf = pd.DataFrame(sgv, columns=SGV._fields)
    fdf = pd.DataFrame(np.array(pop.fixations.array()))
    fdf['ftime'] = pop.fixation_times
    fdf = fdf[fdf.neutral == 0]  # Remove neutral fixations
    fdf.drop(['neutral', 'h'], axis=1, inplace=True)
    fdf = fdf[fdf.ftime > 10 * pop.N]  # Remove fixations prior to shift
    # Remove fixations arising 100 gens after shift
    fdf = fdf[fdf.g <= 10 * pop.N + 100]
    fdf['locus'] = np.array(fdf.pos / 12.0, dtype=np.int32)
    fdf['new'] = [0] * len(fdf.index)
    fdf['standing'] = [0] * len(fdf.index)
    ghat = 2.0 * np.sqrt(2.0 * mu)
    fdf.loc[(fdf.g < 10 * pop.N) & (fdf.s.abs() >= ghat), 'standing'] = 1
    fdf.loc[(fdf.g >= 10 * pop.N) & (fdf.s.abs() >= ghat), 'new'] = 1

    # Change dtype for smaller databases
    fdf.standing = np.array(fdf.standing, dtype=np.int32)
    fdf.new = np.array(fdf.new, dtype=np.int32)
    fdf = fdf.merge(sgvdf, on=['pos', 's', 'g'], how='left',
                    suffixes=('', '_y'))
    fdf['nhaps_at_shift'].fillna(-1, inplace=True)
    fdf.nhaps_at_shift = pd.to_numeric(fdf.nhaps_at_shift,
                                       downcast='integer')
    return fdf


def integrate_genome_scan_with_fixations(pop, fdf, gsdata, repid):
    gsdf = pd.DataFrame(gsdata, columns=GSData._fields)
    gsdf['repid'] = [repid] * len(gsdf.index)

    # We need to SUM within locus/type
    dfg = fdf.groupby(['locus'])['new', 'standing'].sum().reset_index()
    gsdf = gsdf.merge(dfg.ix[:, ['locus', 'standing', 'new']], on=[
        'locus'], suffixes=("", "_y"),
        how='left')
    gsdf[['standing', 'new']] = gsdf[['standing', 'new']].fillna(value=0.0)
    return gsdf


if __name__ == "__main__":
    parser = make_parser()
    args = parser.parse_args(sys.argv[1:])

    # make our named tuple
    qdf = ['generation', 'vg', 'zbar']
    for i in range(args.nloci):
        qdf.append('g' + str(i))
    QData = namedtuple('QData', qdf)

    initial_seed = args.seed
    np.random.seed(initial_seed)
    repseeds = np.random.randint(0, 42000000, args.nreps)
    arglist = [(args, i, j) for i, j in zip(range(len(repseeds)), repseeds)]

    for i in [args.qoutfile, args.goutfile, args.foutfile]:
        if os.path.exists(i):
            os.remove(i)
    mt = args.ncores
    with concurrent.futures.ProcessPoolExecutor(mt) as executor:
        futures = {executor.submit(runsim, i): i for i in arglist}
        for fut in concurrent.futures.as_completed(futures):
            repid, pop, sampler = fut.result()
            fdf = classify_fixations(pop, sampler.sgv, args.mu)
            with sqlite3.connect(args.qoutfile) as conn:
                qdf = pd.DataFrame(sampler.qdata, columns=QData._fields)
                qdf['repid'] = [repid] * len(qdf.index)
                qdf.to_sql('data', conn, if_exists='append', index=False,
                           chunksize=1000)
            with sqlite3.connect(args.goutfile) as conn:
                gsdf = integrate_genome_scan_with_fixations(
                    pop, fdf, sampler.gsdata, repid)
                gsdf.to_sql('data', conn, if_exists='append', index=False,
                            chunksize=1000)
            with sqlite3.connect(args.foutfile) as conn:
                fdf['repid'] = [repid] * len(fdf.index)
                fdf.to_sql('data', conn, if_exists='append',
                           index=False, chunksize=1000)
