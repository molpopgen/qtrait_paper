import pandas as pd
import glob,sys

popstat_files=sorted([i for i in glob.glob('../*stats.h5') if i.find('10regions')==-1])

#fxn to pull simulation params out of file names
def parse_params(fn):
    x=fn.split('_')
    opt=x[3]
    mu=x[5]
    return {'opt':opt,'mu':mu,'sigmu':0.25}

meanSstats=[]

for i in popstat_files:
    pi = parse_params(i)
    x=pd.read_hdf(i)
    xm=x.groupby(['generation','stat']).mean()
    xm.reset_index(inplace=True)
    for k in pi:
        xm[k]=[pi[k]]*len(xm.index)
    meanSstats.append(xm)

out=pd.HDFStore('oneLocusQstatMeans.h5','w',complevel=6,complib='zlib')
out.append('means',pd.DataFrame(pd.concat(meanSstats)))
out.close()

