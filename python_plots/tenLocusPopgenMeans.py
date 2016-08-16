import pandas as pd
import glob

summstat_files=sorted(glob.glob('../*_10regions_popgen.h5'))
popstat_files=sorted(glob.glob('../*_10regions_stats.h5'))

#fxn to pull simulation params out of file names
def parse_params(fn):
    x=fn.split('_')
    opt=x[3]
    mu=x[5]
    return {'opt':opt,'mu':mu,'sigmu':0.25}

meanSstats=[]

for i,j in zip(summstat_files,popstat_files):
    pi = parse_params(i)
    pj = parse_params(j)
    if pi != pj:
        raise RuntimeError("oops: file names must not be sorted!")
    x=pd.read_hdf(i)
    xm=x.groupby(['generation','variable']).mean()
    xm.reset_index(inplace=True)
    for k in pi:
        xm[k]=[pi[k]]*len(xm.index)
    meanSstats.append(xm)

out=pd.HDFStore('tenLocusPopgenMeans.h5','w',complevel=6,complib='zlib')
out.append('means',pd.DataFrame(pd.concat(meanSstats)))
out.close()
