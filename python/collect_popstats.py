import glob
import pandas as pd
import feather

ifiles=glob.glob('*.popstats.h5')

data = [pd.read_hdf(i) for i in ifiles]

for i in range(len(ifiles)):
    fi=ifiles[i].split('_')
#1,3,5
    mu=float("0."+str(fi[5].split('.')[1]))
    data[i]['H2']=[fi[1]]*len(data[i].index)
    data[i]['opt']=[fi[3]]*len(data[i].index)
    data[i]['mu']=[mu]*len(data[i].index)

data=pd.concat(data)

feather.write_dataframe(data,'popstats.feather')
    
