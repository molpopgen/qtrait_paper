##Basic diversity stats
##Outputs a 'tidy' pandas.DataFrame to an H5 file,
##With mean of each summary statistic as function of generation

import libsequence.polytable as polyt
import libsequence.summstats as sstats

try:
    import cPickle as pickle
except ImportError:
    import pickle

import gzip,sys,glob
import pandas as pd

ifile=sys.argv[1]
nbatches=int(sys.argv[2])
ofile=sys.argv[3]

PIK=gzip.open(ifile,"rb")

statNames=['tajd','thetaw','thetapi',"nd1",
           'hprime','nSL','iHS','H12','H1','H2H1']

#Collect all results here

def makeDF(i,j):
    ps=sstats.polySIM(i)
    x=sstats.std_nSLiHS(i)
    g=sstats.garudStats(i)
    values=[ps.tajimasd(),
            ps.thetaw(),
            ps.thetapi(),
            ps.numexternalmutations(),
            ps.hprime(),x[0],x[1],
            g['H12'],g['H1'],g['H2H1']]
    return pd.DataFrame({'gen':[j]*len(values),'stats':statNames,'values':values})

SUMMSTATS=pd.DataFrame()

for i in range(nbatches):
    #There are 3 data objects dumped per sample
    for I in range(3):
        for j in range(64):
            data=pickle.load(PIK)
        #Get "SimData" objects for all neutral sites
            ss=[polyt.simData(k[1]['genotypes'][0]) for k in data]
        #Basic stats
            stats=[makeDF(k,l[0]) for k,l in zip(ss,data)]
            SUMMSTATS=pd.concat([SUMMSTATS,pd.concat([pd.DataFrame(k) for k in stats])])

MEANS=SUMMSTATS.groupby(['gen','stats']).mean()
MEANS.reset_index(inplace=True)

hdf=pd.HDFStore(ofile,'w',complevel=6,complib='zlib')
hdf.open()  
hdf.append('stats',MEANS)
hdf.close()
