import pandas as pd
import glob,sys,gc
import multiprocessing as mp
from qstatmeans import process_file,parse_params,update_totals

popstat_files=sorted([i for i in glob.glob('../*_stats.h5') if i.find('10regions')==-1])
ofilename='oneLocusQstatMeans.h5'
out=pd.HDFStore(ofilename,'w',complevel=6,complib='zlib')
out.close()
P=mp.Pool(10)
P.imap_unordered(process_file,[(ofilename,i,parse_params(i)) for i in popstat_files])
P.close()
P.join()
#for i in popstat_files:
#    pi = parse_params(i)
#    process_file(out,i,pi)
#out.close()

