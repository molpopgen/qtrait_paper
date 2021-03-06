import pandas as pd
import numpy as np
from loci_contributing import *
import multiprocessing

x=pd.read_hdf('tenLocusSweepSummaries.h5')
Fixations = process_fixations(x,True)
OUTFILE='tenLocusPopgenMeansLociWithFixations.h5'
out=pd.HDFStore(OUTFILE,'w',complevel=6,complib='zlib')
out.close()
NFg = Fixations.groupby(['opt','mu'])
params=[]
for n,g in NFg:
    opt=str(n[0])
    if opt=='1.0':
        opt='1'
    mu=str(n[1])
    fn = '../H2_1.0_OPT_' + opt + '_mu_' + mu + '_sigmu0.25_10regions_popgen.h5'
    params.append((g,fn,opt,mu,OUTFILE))

P=multiprocessing.Pool(2)
P.map(process_file,params,chunksize=2)
P.close()
P.join()

