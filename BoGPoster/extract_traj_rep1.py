import pandas as pd
import feather

x=pd.read_hdf('H2_0.2_OPT_0.5_mu_0.001.traj.h5','trajectories')
x=x[x.rep==0]
feather.write_dataframe(x,'H2_0.2_OPT_0.5_mu_0.001.traj.rep1.feather')
