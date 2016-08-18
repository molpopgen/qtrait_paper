
# coding: utf-8

# Extract summary stats of variation at loci that do NOT contribute to adaptation.

# In[1]:

import pandas as pd
import numpy as np


# In[2]:

x=pd.read_hdf('tenLocusSweepSummaries.h5')


# In[3]:

#This is the expected range of loci.
#if any are missing in x for a specific rep, 
#then those loci did not have sweeps from standing varation
#nor any fixations after ('hard sweeps')
loci=set(range(10))


# In[4]:

#build a list of dicts of loci not involved
temp=[]
xg=x.groupby(['opt','mu'])
for n,g in xg:
    reps=g.groupby(['rep'])
    for n2,g2 in reps:
        #Get loci that DID 
        #have fixations for this rep
        l=set(g2.locus)
        nofixations = loci.difference(l)
        for nf in nofixations:
            t = {'opt':n[0],'mu':n[1],'rep':n2,
                'locus':nf}
            temp.append(t)
NoFixations=pd.DataFrame(temp)
            


# Step 2: grouping by optimum/mutation rate, read in the variation data and make some plots

# In[ ]:

out=pd.HDFStore('tenLocusPopgenMeansLociWithoutFixations.h5','w',complevel=6,complib='zlib')
NFg = NoFixations.groupby(['opt','mu'])
for n,g in NFg:
    opt=str(n[0])
    if opt=='1.0':
        opt='1'
    mu=str(n[1])
    fn = '../H2_1.0_OPT_' + opt + '_mu_' + mu + '_sigmu0.25_10regions_popgen.h5'
    d = pd.read_hdf(fn)
    reps=g.groupby(['rep'])
    d.set_index(['rep','locus'],inplace=True)
    temp=pd.DataFrame(g)
    temp.set_index(['rep','locus'],inplace=True)
    result=temp.join([d],how='inner').reset_index()
    means=result.groupby(['generation','variable']).mean().reset_index()
    out.append('nofixations',means)
out.close()
