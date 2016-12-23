import pandas as pd
import glob,sys
from datetime import datetime
from qstatmeans import process_file,parse_params,update_totals
import multiprocessing as mp

summstat_files=sorted(glob.glob('../*_10regions_popgen.h5'))
popstat_files=sorted(glob.glob('../*_10regions_stats.h5'))

for i,j in zip(summstat_files,popstat_files):
    pi = parse_params(i)
    pj = parse_params(j)
    if pi != pj:
        raise RuntimeError("oops: file names must be sorted!")

ofilename='tenLocusQstatMeans.h5'
out=pd.HDFStore(ofilename,'w',complevel=6,complib='zlib')
out.close()

P=mp.Pool(20)
P.imap_unordered(process_file,[(ofilename,i,parse_params(i)) for i in popstat_files])
P.close()
P.join()
