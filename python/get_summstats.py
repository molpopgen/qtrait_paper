##Basic diversity stats
##Outputs a 'tidy' pandas.DataFrame to an H5 file

import libsequence.polytable as polyt
import libsequence.summstats as sstats

try:
    import cPickle as pickle
except ImportError:
    import pickle

import gzip,sys,glob
import pandas as pd

ifiles=sys.argv[1]
nbatches=int(sys.argv[2])
ofile=sys.argv[3]

PIK=gzip.open(ifile,"rb")

statNames=['tajd','thetaw','thetapi',"nd1",
           'hprime']

#Collect all results here

def makeDF(i,j):
    values=[i.tajimasd(),
            i.thetaw(),
            i.thetapi(),
            i.numexternalmutations(),
            i.hprime()]
    print len(values)," ",len(statNames)
    return pd.DataFrame({'gen':[j]*len(values),'stats':statNames,'values':values})

hdf=pd.HDFStore(ofile,'w',complevel=6,complib='zlib')
hdf.open()    

for i in range(nbatches):
    #There are 3 data objects dumped per sample
    for I in range(3):
        for j in range(64):
            data=pickle.load(PIK)
            print len(data)
        #Get "SimData" objects for all neutral sites
            ss=[polyt.simData(k[1]['genotypes'][0]) for k in data]
        #Now, "PolySIM"
            ps=[sstats.polySIM(k) for k in ss]
        #Basic diversity stats
            stats=[makeDF(k,l[0]) for k,l in zip(ps,data)]
            hdf.append('stats',pd.concat([pd.DataFrame(k) for k in stats]))

hdf.close()


