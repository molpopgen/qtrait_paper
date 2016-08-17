import pandas as pd
import glob

x=pd.read_hdf('tenLocusSweepSummaries.h5')

xgParams = x.groupby(['opt','mu'])

LociInvolvedTemp = []
#We need to check something: there are no records in the h5 file for replicates that had no hard/soft sweeps occur!
#We need to find those replicates, and add in a set of 0s for them...
repIDs=range(1024)
for n,g in xgParams:
    dfTempRep=g.groupby(['rep'])
    I=0
    for r,gr in dfTempRep:
        if(r > repIDs[I]): #There is a missing replicate in the h5 file, meaning no fixations (this is rare...)
            v={'opt':n[0],
              'mu':n[1],
              'rep':repIDs[I],
              'nhard':0,
              'nsoft':0,
              'hardLoci':0,
              'softLoci':0,
              'hardAndSoftLoci':0}
            LociInvolvedTemp.append(v)
            I+=1
        softLoci = gr[gr.type=='s'].locus
        hardLoci = gr[gr.type=='h'].locus
        nhardSweeps = len(hardLoci)
        nsoftSweeps = len(softLoci)
        uhardLoci = sorted(hardLoci.unique())
        usoftLoci = sorted(softLoci.unique())
        v={'opt':n[0],
           'mu':n[1],
           'rep':r,
           'nhard':nhardSweeps,
           'nsoft':nsoftSweeps,
           'hardLoci':len(uhardLoci),
           'softLoci':len(usoftLoci),
           'hardAndSoftLoci':len(list(set(uhardLoci)&set(usoftLoci)))}
        LociInvolvedTemp.append(v)
        I+=1

LociInvolved=pd.DataFrame(LociInvolvedTemp)

LociInvolvedGroup = LociInvolved.groupby(['opt','mu']).mean().reset_index()
LociInvolvedGroup.to_csv("NumberLociInvolved.csv",index=False,sep=' ',float_format='%.3f')
