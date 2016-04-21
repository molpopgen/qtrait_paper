import pandas as pd
import numpy as np
import sys
import feather

infile = sys.argv[1]


data1=pd.read_hdf(infile)
#Get allele age for every variant at every generation...
data1['age']=pd.Series(data1.generation-data1.origin+np.array([1]*len(data1.index),dtype=np.float64))
data1.set_index(['origin'],drop=False,inplace=True)
#Select mutations arising N generations b4 optimum shift
data1=data1[data1.origin>9000]
#Group data
g=data1.groupby(['origin','pos','esize','rep'])

##Collect mutations that were lost
losses = g.filter(lambda x:x.freq.max()<1.0)
##Write these to a file
ofile=infile.replace('.traj.h5','.losses.h5')
losses.to_hdf(ofile,'losses',complevel=6,complib='zlib')
##Get mean values of everything per generation
losses=losses.groupby('generation').mean()
ofile=infile.replace('.traj.h5','.losses_means.h5')
losses.to_hdf(ofile,'losses_means',complevel=6,complib='zlib')

del losses

#Not do fixations
fixations = g.filter(lambda x:x.freq.max()==1.0)
ofile=infile.replace('.traj.h5','.fixations.h5')
fixations.to_hdf(ofile,'fixations',complevel=6,complib='zlib')
fixations=fixations.groupby('generation').mean()
ofile=infile.replace('.traj.h5','.fixations_means.h5')
fixations.to_hdf(ofile,'fixations_means',complevel=6,complib='zlib')
