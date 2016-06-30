from __future__ import print_function
import fwdpy as fp
import fwdpy.qtrait_mloc as qtm
import libsequence.polytable as polyt
import libsequence.extensions as ext
import KRT_qtrait_paper.summstatsParallel as sstats
import numpy as np
import pandas as pd
import getopt
import sys

def usage():
    print ("something went wrong")

valid_trait_models=['additive','mult']
trait_models = {'additive':qtm.MlocusAdditiveTrait(),
                'mult':qtm.MlocusMultTrait()}

valid_sampler_names=['stats','freq','popgen']

def get_sampler(samplerString,length,optimum,nsam,rng):
    if samplerString == 'stats':
        return fp.QtraitStatsSampler(length,optimum)
    elif samplerString == 'freq':
        return fp.FreqSampler(length)
    elif samplerString == 'popgen':
        if nsam is None:
            print("sample size cannot be none when sampler is ",sampler_string)
        return fp.PopSampler(length,nsam,rng)
    else:
        print ("invalid sampler name")
        usage()
        sys.exit(0)

def get_trait_model(traitString):
   return trait_models[traitString]

def write_output(sampler,output,nloci,REPID):
    if isinstance(sampler,fp.FreqSampler):
        data=[pd.DataFrame(i) for i in fp.tidy_trajectories(sampler.get())]
        for df in data:
            df['rep']=[REPID]*len(df.index)
            REPID+=1
        output.append('trajectories',pd.concat(data))
    elif isinstance(sampler,fp.QtraitStatsSampler):
        data=[pd.DataFrame(i) for i in sampler.get()]
        for df in data:
            df['rep']=[REPID]*len(df.index)
            REPID+=1
        output.append('stats',pd.concat(data))
    elif isinstance(sampler,fp.PopSampler):
        data=[pd.DataFrame(i) for i in sampler.get()]
        for REP in data:
            for LOCUS in range(nloci):
                gen=[i[0] for i in REP]
                locus=[i[1][LOCUS]['genotypes'][0] for i in REP]
                v=ext.SimDataVec(locus)
                stats = sstats.getSummStatsParallel(v)
                DF=[pd.DataFrame(i.items(),columns=['stat','value']) for i in stats]
                for i in range(len(DF)):
                    DF[i]['gen']=[gen[i]]*len(DF[i].index)
                    DF[i]['locus']=[LOCUS]*len(DF[i].index)
                    DF[i]['replicate']=[REPID]*len(DF[i].index)
                    rv.append(pd.concat(DF))
            REPID+=1
    else:
        raise RuntimeError("uh oh: sampler type not recognized for output.  We shouldn't have gotten this far!")
    return REPID

def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:],"m:e:H:S:O:N:t:s:F:r:n:d:",
                                   ["theta=","rho=","trait=","sampler=","nsam=","cores=","batches=",
                                    "nloci="])
    except getopt.GetoptError as err:
        # print help information and exit:
        print(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
    #set up default params
    N=1000   # pop size
    t=None     # 0.1N
    e = 0.25 # s.d. of effect sizes
    S = 1    # V(S)
    H = 1.0 # desired b-sense H^2
    m = None # Mutation rate (per gamete, per generation) to alleles affecting trait value
    r = 0.5 #Recombination rate between locis
    Opt = 0.0  # Value of optimum after 10N gens
    ofile=None
    seed = 0
    NLOCI=10
    theta = 100 # mutational size of region where neutral mutations occur
    rho = 100   # recombination size of region where neutral mutations occur
    traitString="additive"
    samplerString=None
    ncores=64
    nbatches=16
    dominance=1.0
    nsam = None
    for o,a in opts:
        if o == '-m':
            m = float(a)
        elif o == '-e':
            e = float(a)
        elif o == '-H':
            H = float(a)
        elif o == '-S':
            S = float(a)
        elif o == '-O':
            Opt = float(a)
        elif o == '-N':
            N=int(a)
        elif o == '-t':
            t = int(a)
        elif o == '-s':
            seed = int(a)
        elif o == '-F':
            ofile = a
        elif o == '-r':
            r = float(a)
        elif o == '--nsam':
            ssize = int(a)
        elif o == '--theta':
            theta=float(a)
        elif o == '--rho':
            rho=float(a)
        elif o == '--trait':
            traitString = a
            if a not in valid_trait_models:
                print ("Error: invalid trait model")
                usage()
                sys.exit(0)
        elif o == '--sampler':
            samplerString = a
            if a not in valid_sampler_names:
                print ("Error: invalid sampler name")
                usage()
                sys.exit(0)
        elif o == '--cores':
            ncores=int(a)
            if ncores < 1:
                print ("--ncores must be > 0")
                usage()
                sys.exit(0)
        elif o == '--batches':
            nbatches=int(a)
            if nbatches < 1:
                print ("--nbatches myst be > 0")
                usage()
                sys.exit(0)
        elif o == '--nloci':
            NLOCI=int(a)
            if NLOCI<2:
                rpint ("--nloci must be > 1")
                usage()
                sys.exit(0)

    if samplerString is None:
        print ("Error: sampler must be defined")
        usage()
        sys.exit(0)

    if t is None:
        print ("Error: sampling interval must be defined")
        usage()
        sys.exit(0)
    if m is None:
        usage()
        sys.exit(2)
    if ofile is None:
        usage()
        sys.exit(2)

    #Can start working now:
    REP=0
    out=pd.HDFStore(ofile,"w",complevel=6,complib='zlib')

    little_r_per_locus = rho/(4.0*float(N))
    mu_n_region=theta/(4.0*float(N))

    rnge=fp.GSLrng(seed)
    rngs=fp.GSLrng(seed)
    nlist=np.array([N]*(10*N),dtype=np.uint32)
    fitness = get_trait_model(traitString)
    sregions=[fp.GaussianS(0,1,1,e,dominance)]*NLOCI
    for BATCH in range(nbatches):
        x = fp.MlocusPopVec(ncores,N,NLOCI)

        sampler=get_sampler(samplerString,len(x),Opt,nsam,rngs)
        qtm.evolve_qtraits_mloc_sample_fitness(rnge,x,sampler,fitness,nlist,
                                               [mu_n_region]*NLOCI,
                                               [m/float(NLOCI)]*NLOCI,
                                               sregions,
                                               [little_r_per_locus]*NLOCI,
                                               [r]*(NLOCI-1),#loci unlinked
                                               sample=t,VS=S)
        REP=write_output(sampler,out,NLOCI,REP)
        sampler=get_sampler(samplerString,len(x),Opt,nsam,rngs)
        qtm.evolve_qtraits_mloc_sample_fitness(rnge,x,sampler,fitness,nlist,
                                               [mu_n_region]*NLOCI,
                                               [m/float(NLOCI)]*NLOCI,
                                               sregions,
                                               [little_r_per_locus]*NLOCI,
                                               [r]*(NLOCI-1),#loci unlinked
                                               sample=t,VS=S,optimum=Opt)
        REP=write_output(sampler,out,NLOCI,REP)

    out.close()

if __name__ == "__main__":
    main()
