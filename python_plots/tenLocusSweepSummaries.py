import glob
import pandas as pd
from common_functions import *
files=glob.glob('../*10regions_popgen_fixations.h5')
dflist=[]
for f in files:
    fi=parse_params(f)
    x=pd.read_hdf(f)
    x=x[x.neutral==False]

    #separate sweeps from standing variation
    #from sweeps from new mutations
    hard = x[x.g >= 100000].reset_index()
    soft = x[(x.g < 100000) & (x.ftime >= 100000)].reset_index()
    hard['type']=['h']*len(hard.index)
    soft['type']=['s']*len(soft.index)
    for i in fi:
        hard[i]=[fi[i]]*len(hard.index)
        soft[i]=[fi[i]]*len(soft.index)
    dflist.append(pd.DataFrame(pd.concat([hard,soft])))

df=pd.DataFrame(pd.concat(dflist))

out=pd.HDFStore('tenLocusSweepSummaries.h5','w')
out.append('sweeps',df)
out.close()
