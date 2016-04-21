/usr/bin/time -f "%e %M" -o reduce1.txt python reduce_trajectory_data.py ../H2_0.2_OPT_1_mu_0.001.traj.h5 > reduce1.stdout 2> reduce1.stderr & 
/usr/bin/time -f "%e %M" -o reduce2.txt python reduce_trajectory_data.py ../H2_0.2_OPT_0.5_mu_0.001.traj.h5  > reduce2.stdout 2> reduce2.stderr &
