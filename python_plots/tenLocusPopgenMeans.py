import pandas as pd
import glob
from common_functions import *
import multiprocessing as mp

summstat_files=sorted(glob.glob('../*_10regions_popgen.h5'))
popstat_files=sorted(glob.glob('../*_10regions_stats.h5'))
def update_totals(args):
    filei=args[0]
    Totals=args[1]
    START=args[2]
    STOP=args[3]
    print filei,START,STOP
    x=pd.read_hdf(filei,where='rep>=START&rep<STOP')
    x.reset_index(inplace=True)
    x.drop_duplicates(inplace=True)
    xm=x.groupby(['generation','variable']).sum().reset_index()
    if Totals is None:
        Totals=xm
    else:
        temp = Totals.append(xm)
        tempg=temp.groupby(['generation','variable']).sum().reset_index()
        Totals=temp
    return Totals

def process_file(ifile,params,out):
    print "processing ",ifile
    START=0
    STOP=64
    Totals=None
    while START < 1024:
        P=mp.Pool(1)
        Totals=P.map(update_totals,[(ifile,Totals,START,STOP)])[0]
        P.close()
        #x=pd.read_hdf(ifile,where='rep >= START & rep < STOP')
        #x.reset_index(inplace=True)
        #x.drop_duplicates(inplace=True)
        #xm=x.groupby(['generation','variable']).mean()
        #xm.reset_index(inplace=True)
        START+=64
        STOP+=64
    for k in params:
        Totals[k]=[params[k]]*len(Totals.index)
    Totals.value/=float(10*1024/64)
    del Totals['rep']
    del Totals['locus']
    out.append('means',Totals)

out=pd.HDFStore('tenLocusPopgenMeans.h5','w',complevel=6,complib='zlib')

for i,j in zip(summstat_files,popstat_files):
    pi = parse_params(i)
    pj = parse_params(j)
    if pi != pj:
        raise RuntimeError("oops: file names must not be sorted!")
    print "reading ",i
    process_file(i,pi,out)
#    START=0
#    STOP=64
#    while START < 1024:
#        x=pd.read_hdf(i,where='rep >= START & rep < STOP')
#        x.reset_index(inplace=True)
#        x.drop_duplicates(inplace=True)
#        xm=x.groupby(['generation','variable']).mean()
#        xm.reset_index(inplace=True)
#        for k in pi:
#            xm[k]=[pi[k]]*len(xm.index)
#        del xm['rep']
#        del xm['locus']
#        running_means.append(xm)
#        START+=64
#        STOP+=64
#df=pd.concat(running_means)
#xm=df.groupby(['generation','variable']).mean().reset_index()
#out.append('means',xm)  
out.close()
