import pandas as pd
import glob,sys
from datetime import datetime
from qstatmeans import process_file,parse_params,update_totals
summstat_files=sorted(glob.glob('../*_10regions_popgen.h5'))
popstat_files=sorted(glob.glob('../*_10regions_stats.h5'))

out=pd.HDFStore('tenLocusQstatMeans.h5','w',complevel=6,complib='zlib')
for i,j in zip(summstat_files,popstat_files):
    pi = parse_params(i)
    pj = parse_params(j)
    if pi != pj:
        raise RuntimeError("oops: file names must be sorted!")
    process_file(out,j,pj)
out.close()
