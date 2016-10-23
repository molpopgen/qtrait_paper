#Simulate the Tennessen model with the following params:
#1. 10 loci separated by 50cM
#2. Each locus is 11*theta and 11*rho long, where theta & rho are per-locus input params
#3. The middle of each locus is where causative mutations happen at rate mu with 
#   and sigmamu=0.25 and VS=1
#4. A sample of n=6e3 is taken at the end and summary stats are calcuated for 11 window
#5. The samples are written to files for latter analysis with SweeD and SDS
#6. The change in the optimum trait value occurs "out of Africa".  The optimum is 0 prior to then.

import fwdpy as fp
import fwdpy.qtrait_mloc as qtm
import fwdpy.demography as fpd
import numpy as np
import pandas as pd
import sys
import datetime
import getopt
import warnings
import libsequence.parallel as lsp
import libsequence.polytable as pt
import libsequence.summstats as sstats
import libsequence.windows as windows

def get_nlist1():
    """
    Generates a numpy array of the canges in N over time

    There are 5 epochs, with t=0 being the present.
    
    E1: Ne= 7,310  from t=start(8N?) to t = - 5920 generation (Ancestral sizes until 5920 generations ago)
    
    E2: Ne =14,474 from t = -5920 to t = -2040 (Ancient growth at 5920 g ago)
    
    E3: Ne =1,861 from t= -2040 to t= -920 (OOA, first bottle neck 2040 g ago)
    
    E4: Ne = 1,032 to Ne = 9,300 during t = -920 to t = -205 ( second bottle neck and onset of 715 g of exponential growth at rate 0.31% per gen )  
    
    E5: Ne = 9,300 to Ne = 512,000 during t = -205 to t = -0 ( 205 g of exponential growth at rate 1.95% per gen )  
    """
    n=[7310]*(10*7310) #E1: evolve ancestral size to mutation/selection/drift equilibrium
    n.extend([14474]*(5920-2040)) #E2
    #n.extend([1861]*(2040-920)) #E3
    #n.extend(fpd.exponential_size_change(1032,9300,920-205)) #E4
    #n.extend(fpd.exponential_size_change(9300,51200,205)) #E5
    return n

def get_nlist2():
    """
    Generates a numpy array of the canges in N over time

    There are 5 epochs, with t=0 being the present.
    
    E1: Ne= 7,310  from t=start(8N?) to t = - 5920 generation (Ancestral sizes until 5920 generations ago)
    
    E2: Ne =14,474 from t = -5920 to t = -2040 (Ancient growth at 5920 g ago)
    
    E3: Ne =1,861 from t= -2040 to t= -920 (OOA, first bottle neck 2040 g ago)
    
    E4: Ne = 1,032 to Ne = 9,300 during t = -920 to t = -205 ( second bottle neck and onset of 715 g of exponential growth at rate 0.31% per gen )  
    
    E5: Ne = 9,300 to Ne = 512,000 during t = -205 to t = -0 ( 205 g of exponential growth at rate 1.95% per gen )  
    """
    n=[1861]*(2040-920) #E3
    n.extend(fpd.exponential_size_change(1032,9300,920-205)) #E4
    n.extend(fpd.exponential_size_change(9300,51200,205)) #E5
    return n

def main():
    #This is the number of loci where causative variants occur.
    #Each locus will be in the "middle" of a region, so there
    #will be 10x more DNA on either side
    NLOCI=10
    NREPS=1
    NCORES=1
    SEED=None
    SIGMU=0.25
    TBB=-1
    mu=1e-3
    optimum = None
    #Mutation/rec params
    theta_neutral=100.0
    rho_per_locus=100.0
    statfile = None
    try:
        opts,args=getopt.getopt(sys.argv[1:],"",['theta=','rho=','seed=','ncores=','nreps=','tbb=','mu=','statfile=','opt='])
    except getopt.GetoptError as err:
        print err
        sys.exit(1)
    
    for o,a in opts:
        if o == '--theta':
            theta_neutral=float(a)
        elif o == '--rho':
            rho_per_locus=float(a)
        elif o == '--seed':
            SEED=int(a)
        elif o == '--ncores':
            NCORES=int(a)
        elif o == '--nreps':
            NREPS=int(a)
        elif o == '--tbb':
            TBB=int(a)
        elif o == '--mu':
            if(a <= 0):
                raise RuntimeError("mutation rate must be > 0")
            mu=float(a)
        elif o == '--statfile':
            statfile=a
        elif o == '--opt':
            optimum=float(a)

    if optimum is None:
        raise RuntimeError("optimum trait value not specified, use --opt <optimum>")

    if SEED is None:
        warnings.warn("Setting seed to 0, use --seed <seedval> to change")
        SEED=0
    if statfile is None:
        warnings.warn("No output file for summary stats specified.  Use --statfile <filename.h5> to change.")

    nstub = "neutral.mu"+str(mu)+".opt"+str(optimum)
    sstub = "selected.mu"+str(mu)+".opt"+str(optimum)
    rnge=fp.GSLrng(SEED)

    START=2500
    NANC=7310
    locus_boundaries=[(float(i+i*11),float(i+i*11+11)) for i in range(NLOCI)]
    nregions=[fp.Region(j[0],j[1],theta_neutral/(4.*float(NANC)),coupled=True) for i,j in zip(range(NLOCI),locus_boundaries)]
    recregions=[fp.Region(j[0],j[1],rho_per_locus/(4.*float(NANC)),coupled=True) for i,j in zip(range(NLOCI),locus_boundaries)]
    sregions=[fp.GaussianS(j[0]+5.,j[0]+6.,mu,SIGMU,coupled=False) for i,j in zip(range(NLOCI),locus_boundaries)]
    f=qtm.MlocusAdditiveTrait()
    sched = lsp.scheduler_init(TBB)
    statnames=['S','tajimasd','pi','thetaw','hprime','H1','H12','H2H1','nSL','iHS']
    H5out = None
    if statfile is not None:
        H5out = pd.HDFStore(statfile,'w')
    repid=0
    for batch in range(NREPS):
        nlist=np.array(get_nlist1(),dtype=np.uint32)
        #nlist=np.array([NANC]*100,dtype=np.uint32)
        pops = fp.MlocusPopVec(NCORES,nlist[0],NLOCI)
        sampler=fp.NothingSampler(len(pops))
        d=datetime.datetime.now()
        print d.now()
        qtm.evolve_qtraits_mloc_regions_sample_fitness(rnge,pops,sampler,f,
                nlist[0:],
                nregions,sregions,recregions,
                [0.5]*(NLOCI-1),
                0,
                0,0.)
        d=datetime.datetime.now()
        print d.now()
        nlist=np.array(get_nlist2(),dtype=np.uint32)
        #nlist=np.array([NANC]*100,dtype=np.uint32)
        qtm.evolve_qtraits_mloc_regions_sample_fitness(rnge,pops,sampler,f,
                nlist[0:],
                nregions,sregions,recregions,
                [0.5]*(NLOCI-1),
                0,optimum)
        d=datetime.datetime.now()
        print d.now()
        neutralFile = nstub + '.batch' + str(batch)
        selectedFile = sstub + '.batch' + str(batch)
        sampler=fp.PopSampler(len(pops),6000,rnge,False,neutralFile,selectedFile,recordSamples=True,boundaries=locus_boundaries)
        fp.apply_sampler(pops,sampler)
        d=datetime.datetime.now()
        print d.now()
        if statfile is not None:
            dflist=[]
            #get data from sampler.  A list of generators is returned
            data = sampler.get()
            #get summary stats in sliding windows based on neutral diversity
            for di in data:
                simDataList=[pt.simData(k[0][0]) for k in di.yield_data()]
                w=[windows.Windows(i,window_size=1.,step_len=1.,starting_pos=j[0],ending_pos=j[1]) for i,j in zip(simDataList,locus_boundaries)]
                polySIMlist=[]
                hapstats=[] #nSL, etc.
                for i in w:
                    polySIMlist.extend([sstats.polySIM(j) for j in i])
                    hapstats.extend([(sstats.garudStats(j),sstats.std_nSLiHS(j,0.05,0.1)) for j in i])
                stats=[(i.numpoly(),i.tajimasd(),i.thetapi(),i.thetaw(),i.hprime()) for i in polySIMlist]
                combinedStats=[i+(j[0]['H1'],j[0]['H12'],j[0]['H2H1'],j[1][0],j[1][1]) for i,j in zip(stats,hapstats)]
                window=0
                for cs in combinedStats:
                    tempDF=pd.DataFrame({'stat':statnames,'value':list(cs),
                            'window':[window]*len(statnames),'rep':[repid]*len(statnames)})
                    window+=1
                    dflist.append(tempDF)
                repid+=1 
            H5out.append('stats',pd.concat(dflist))
    if statfile is not None:
        H5out.close()

if __name__ == "__main__":
    main()



