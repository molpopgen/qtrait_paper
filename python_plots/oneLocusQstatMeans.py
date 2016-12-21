import pandas as pd
import glob,sys,gc
from qstatmeans import process_file,parse_params,update_totals

popstat_files=sorted([i for i in glob.glob('../*_stats.h5') if i.find('10regions')==-1])

out=pd.HDFStore('oneLocusQstatMeans.h5','w',complevel=6,complib='zlib')
for i in popstat_files:
    pi = parse_params(i)
    process_file(out,i,pi)
out.close()

