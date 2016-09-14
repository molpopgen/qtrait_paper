import glob
import pandas as pd
import numpy as np
from common_functions import *
files=glob.glob('../*10regions_popgen_fixations.h5')
out=pd.HDFStore('tenLocusSweepSummaries.h5','w')
#dflist=[]
for f in files:
    fi=parse_params(f)
    x=pd.read_hdf(f)
    x.reset_index(inplace=True)
    x=x[x.neutral==False]
    #separate sweeps from standing variation
    #from sweeps from new mutations
    hard = x[x.g >= 50000].reset_index()
    soft = x[(x.g < 50000) & (x.ftime >= 50000)].reset_index()
    hard['type']=['h']*len(hard.index)
    soft['type']=['s']*len(soft.index)
    hard['locus']=np.floor(hard.pos)
    soft['locus']=np.floor(soft.pos)
    for i in fi:
        hard[i]=[float(fi[i])]*len(hard.index)
        soft[i]=[float(fi[i])]*len(soft.index)
    out.append('sweeps',pd.concat([hard,soft]))
    #dflist.append(pd.DataFrame(pd.concat([hard,soft])))

#df=pd.DataFrame(pd.concat(dflist))
#out.append('sweeps',df)
out.close()
