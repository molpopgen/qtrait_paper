import pandas as pd

x=pd.read_hdf("../validation_popstats.h5")
means=x.groupby(['mu','r','sige','sigmu']).mean()
means.reset_index(inplace=True)
means=means.rename(columns={'VG':'meanVG'})
stds=x.groupby(['mu','r','sige','sigmu']).std()
stds.reset_index(inplace=True)
means['stdVG']=stds.VG
means['EVG']=4.0*means.mu
del means['rep']

o=pd.HDFStore('validation_summaries.h5',complevel=6,complib='zlib')
o.append('summaries',means)
o.close()
