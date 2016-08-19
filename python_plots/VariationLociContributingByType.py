#Breaks up the summaries in terms of hard vs soft sweeps

import pandas as pd
import numpy as np
from loci_contributing import *

x=pd.read_hdf('tenLocusSweepSummaries.h5')

out=pd.HDFStore('tenLocusPopgenMeansSoftSweepsOnly.h5','w',complevel=6,complib='zlib')
counts = count_soft_hard(x)
counts = counts[(counts.soft>0)&(counts.hard==0)]
cg=counts.groupby(['opt','mu'])
for n,g in cg:
    opt=str(n[0])
    if opt=='1.0':
        opt='1'
    mu=str(n[1])
    fn = '../H2_1.0_OPT_' + opt + '_mu_' + mu + '_sigmu0.25_10regions_popgen.h5'
    d = pd.read_hdf(fn)
    d.set_index(['rep','locus'],inplace=True)
    temp=pd.DataFrame(g)
    temp.set_index(['rep','locus'],inplace=True)
    result=temp.join([d],how='inner').reset_index()
    means=result.groupby(['generation','variable']).mean().reset_index()
    means['opt']=[opt]*len(means.index)
    means['mu']=[mu]*len(means.index)
    means.drop(['rep','locus','soft','hard'],axis=1,inplace=True)
    out.append('fixations',means)
out.close()
        
