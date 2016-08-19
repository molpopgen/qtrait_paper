
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

out=pd.HDFStore('tenLocusPopgenLociWithoutFixations.h5','w',complevel=6,complib='zlib')
NFg = NoFixations.groupby(['opt','mu'])
for n,g in NFg:
    opt=str(n[0])
    if opt=='1.0':
        opt='1'
    mu=str(n[1])
    fn = '../H2_1.0_OPT_' + opt + '_mu_' + mu + '_sigmu0.25_10regions_popgen.h5'
    d = pd.read_hdf(fn)
    reps=g.groupby(['rep'])
    for n2,g2 in reps:
        #We select on generation to reduce output file size.
        #This threshold is arbitrary, but this is still 
        #during the "burn-in" period
        temp=pd.DataFrame(d[(d.generation > 80000) & (d.rep==n2)&(d.locus.isin(g2.locus))])
        temp['opt']=[n[0]]*len(temp.index)
        temp['mu']=[n[1]]*len(temp.index)
        out.append('nofixations',temp)
out.close()
        


# In[ ]:

x=pd.read_hdf('tenLocusPopgenLociWithoutFixations.h5')
xg=x.groupby(['opt','mu']).mean().reset_index()
xg.to_hdf('tenLocusPopgenMeansLociWithoutFixations.h5','means',complevel=6,complib='zlib')

