# Read in the popstats files,
# calculate mean of each stat per
# generation, collate, and output 
# for plotting in R

import glob
import pandas as pd
import feather

ifiles=glob.glob('*.popstats.h5')

means=[]

for i in range(len(ifiles)):
    data=pd.read_hdf(ifiles[i])
    fi=ifiles[i].split('_')
#1,3,5
    data=data.groupby(['generation','stat']).mean()
    data.reset_index(inplace=True)
    mu=float("0."+str(fi[5].split('.')[1]))
    data['H2']=[fi[1]]*len(data.index)
    data['opt']=[fi[3]]*len(data.index)
    data['mu']=[mu]*len(data.index)
    means.append(data)

means=pd.concat(means)
feather.write_dataframe(means,'popstats.feather')


    
