
# coding: utf-8

# Extract summary stats of variation at loci that do NOT contribute to adaptation.

# In[1]:

import pandas as pd
import numpy as np
from loci_contributing import *

# In[2]:

x=pd.read_hdf('tenLocusSweepSummaries.h5')


# In[3]:


Fixations = process_fixations(x,True)

# Step 2: grouping by optimum/mutation rate, read in the variation data and make some plots

# In[ ]:

out=pd.HDFStore('tenLocusPopgenMeansLociWithFixations.h5','w',complevel=6,complib='zlib')
NFg = Fixations.groupby(['opt','mu'])
for n,g in NFg:
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
    out.append('fixations',means)
out.close()
        
