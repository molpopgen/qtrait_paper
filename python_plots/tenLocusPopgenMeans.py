import pandas as pd
import glob
from common_functions import *
summstat_files=sorted(glob.glob('../*_10regions_popgen.h5'))
popstat_files=sorted(glob.glob('../*_10regions_stats.h5'))


#meanSstats=[]
out=pd.HDFStore('tenLocusPopgenMeans.h5','w',complevel=6,complib='zlib')
#out.append('means',pd.DataFrame(pd.concat(meanSstats)))
running_means=[]
for i,j in zip(summstat_files,popstat_files):
    pi = parse_params(i)
    pj = parse_params(j)
    if pi != pj:
        raise RuntimeError("oops: file names must not be sorted!")
    print "reading ",i
    START=0
    STOP=64
    while START < 1024:
        x=pd.read_hdf(i,where='rep >= START & rep < STOP')
        x.reset_index(inplace=True)
        x.drop_duplicates(inplace=True)
        xm=x.groupby(['generation','variable']).mean()
        xm.reset_index(inplace=True)
        for k in pi:
            xm[k]=[pi[k]]*len(xm.index)
        del xm['rep']
        del xm['locus']
        running_means.append(xm)
        START+=64
        STOP+=64
df=pd.concat(running_means)
xm=df.groupby(['generation','variable']).mean().reset_index()
out.append('means',xm)  
out.close()
