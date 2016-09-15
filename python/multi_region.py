from __future__ import print_function
import pyximport
pyximport.install()
import summstatsParallel as plugin
import mlocAges
import PopstatsLocus
import fwdpy as fp
import fwdpy.qtrait_mloc as qtm
import numpy as np
import pandas as pd
import getopt
import sys,math

def usage():
    print ("something went wrong")

valid_trait_models=['additive','mult']
trait_models = {'additive':qtm.MlocusAdditiveTrait(),
                'mult':qtm.MlocusMultTrait()}

valid_sampler_names=['lstats','stats','ages','popgen','ms']

def get_sampler(samplerString,length,optimum,nsam,rng,nstub,sstub):
    if samplerString == 'lstats':
        return PopstatsLocus.PopstatsLocus(length)
    elif samplerString == 'stats':
        return fp.QtraitStatsSampler(length,optimum)
    elif samplerString == 'ages':
        return mlocAges.MlocusAgeSampler(length) 
    elif samplerString == 'popgen':
        if nsam is None:
            print("sample size cannot be none when sampler is ",samplerString)
        return plugin.MlocusSummStatsSampler(length,nsam,rng)
    elif samplerString == 'ms':
        if nsam is None:
            print("sample size cannot be none when sampler is ",samplerString)
#        if nstub is None or sstub is None:
#            raise RuntimeError("file name prefixes cannot be None")
        return fp.PopSampler(length,nsam,rng,False,nstub,sstub,[(float(i),float(i)+1.) for i in
            range(10)],recordSamples=False)
    else:
        print ("invalid sampler name")
        usage()
        sys.exit(0)

def get_trait_model(traitString):
   return trait_models[traitString]

def write_output(sampler,output,nloci,REPID):
    if isinstance(sampler,PopstatsLocus.PopstatsLocus):
        data=[pd.DataFrame(i) for i in sampler.get()]
        for df in data:
            df['rep']=[REPID]*len(df.index)
            REPID+=1
        output.append('lstats',pd.concat(data))
    elif isinstance(sampler,mlocAges.MlocusAgeSampler):
        data=[pd.DataFrame(i) for i in sampler.get()]
        for df in data:
            df['rep']=[REPID]*len(df.index)
            REPID+=1
        output.append('ages',pd.concat(data))
    elif isinstance(sampler,fp.QtraitStatsSampler):
        data=[pd.DataFrame(i) for i in sampler.get()]
        for df in data:
            df['rep']=[REPID]*len(df.index)
            REPID+=1
        output.append('stats',pd.concat(data))
    elif isinstance(sampler,plugin.MlocusSummStatsSampler):
        data=[pd.DataFrame(i) for i in sampler.get()]
        for i in data:
            i['rep']=[REPID]*len(i.index)
            i = pd.melt(i,value_vars=['H1','H12','H2H1','hprime','iHS','nSL','nd1','tajd','thetapi','thetaw'],
                        id_vars=['rep','locus','generation'])
            output.append('summstats',i)
            REPID+=1
    elif isinstance(sampler,fp.PopSampler):
        data=sampler.get()
        for i in data:
            for j in i.yield_data():
                print (j[1])
        #for j in data[0]:
        #    print (len(j))
        #    print (j[1])
        #    t = [ji for ji in j[1] if math.isnan(ji['p'][0]) is False]
        #    print (t)
    else:
        raise RuntimeError("uh oh: sampler type not recognized for output.  We shouldn't have gotten this far!")

def write_fixations(pops,fixationsFileName,repid):
    x=fp.view_fixations(pops)
    hdf=pd.HDFStore(fixationsFileName,'a',complevel=6,complib='zlib'))
    df=[pd.DataFrame(i) for i in x]
    for i in df:
        i['rep']=[repid]*len(i.index)
        repid+=1
        hdf.append('fixations',i)
    hdf.close()

def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:],"m:e:H:S:O:N:t:s:F:r:n:d:",
                                   ["theta=","rho=","trait=","sampler=","nsam=","cores=","batches=",
                                       "nstub=","sstub=",
                                    "nloci=","fixations=","t2=","g2="])#,"neutral="])
    except getopt.GetoptError as err:
        # print help information and exit:
        print(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
    #set up default params
    N=1000   # pop size
    t=None     # 0.1N
    t2=None
    e = 0.25 # s.d. of effect sizes
    S = 1    # V(S)
    H = 1.0 # desired b-sense H^2
    m = None # Mutation rate (per gamete, per generation) to alleles affecting trait value
    r = 0.5 #Recombination rate between locis
    Opt = 0.0  # Value of optimum after 10N gens
    ofile=None
    seed = 0
    NLOCI=10  #The total number of loci to simulate
    #NeutralLocus = 5 #The index of the locus where no causal mutations will occur
    theta = 100 # mutational size of region where neutral mutations occur
    rho = 100   # recombination size of region where neutral mutations occur
    traitString="additive"
    samplerString=None
    ncores=64
    nbatches=16
    dominance=1.0
    ssize = None
    fixationsFileName = None
    G2 = None
    nstub=None
    sstub=None
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
            t2=t
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
        elif o == "--fixations":
            fixationsFileName=a  
        elif o == "--t2":
            t2=int(a)
	elif o == "--g2":
	    G2=int(a)
        elif o == "--nstub":
            nstub=a
        elif o == '--sstub':
            sstub=a
        #elif o == "--neutral":
        #    if a != "None":
        #        NeutralLocus=int(a)
        #    else:
        #        NeutralLocus=a
    if samplerString is None:
        print ("Error: sampler must be defined")
        usage()
        sys.exit(0)

    if t is None:
        print ("Error: sampling interval must be defined")
        usage()
        sys.exit(0)
    if t2 is None:
        t2=t
    if m is None:
        usage()
        sys.exit(2)
    if ofile is None:
        usage()
        sys.exit(2)
    if G2 is None:
	G2=10*N

    #Can start working now:
    REP=0
    out=pd.HDFStore(ofile,"w",complevel=6,complib='zlib')

    if fixationsFileName is not None:
        tfile=pd.HDFStore(fixationsFileName,'w')
        tfile.close()
        
    little_r_per_locus = rho/(4.0*float(N))
    mu_n_region=theta/(4.0*float(N))

    rnge=fp.GSLrng(seed)
    rngs=fp.GSLrng(seed)
    nlist=np.array([N]*(10*N),dtype=np.uint32)
    fitness = get_trait_model(traitString)
    sregions=[fp.GaussianS(0,1,1,e,dominance)]*NLOCI

    neutral_mut_rates = [mu_n_region]*NLOCI
    causal_mut_rates = [m/float(NLOCI)]*NLOCI
    #if NeutralLocus is not None:
    #    causal_mut_rates[NeutralLocus]=0.0
    for BATCH in range(nbatches):
        x = fp.MlocusPopVec(ncores,N,NLOCI)
        nstub_t = None
        sstub_t = None
        if nstub is not None:
            nstub_t = nstub + b'.batch' + str(BATCH)
        if sstub is not None:
            sstub_t = sstub + b'.batch' + str(BATCH)
        sampler=get_sampler(samplerString,len(x),Opt,ssize,rngs,nstub_t,sstub_t)
        qtm.evolve_qtraits_mloc_sample_fitness(rnge,x,sampler,fitness,nlist,
                neutral_mut_rates,
                causal_mut_rates,
                sregions,
                [little_r_per_locus]*NLOCI,
                [r]*(NLOCI-1),#loci unlinked
                sample=t,VS=S)
        sys.exit(1)
        write_output(sampler,out,NLOCI,REP)
        sampler=get_sampler(samplerString,len(x),Opt,ssize,rngs,nstub_t,sstub_t)
        qtm.evolve_qtraits_mloc_sample_fitness(rnge,x,sampler,fitness,nlist[:G2],
                neutral_mut_rates,
                causal_mut_rates,
                sregions,
                [little_r_per_locus]*NLOCI,
                [r]*(NLOCI-1),#loci unlinked
                sample=t2,VS=S,optimum=Opt)
        write_output(sampler,out,NLOCI,REP)
        if fixationsFileName is not None:
            write_fixations(x,fixationsFileName,REP)
        REP+=len(x)

    out.close()

if __name__ == "__main__":
    main()
