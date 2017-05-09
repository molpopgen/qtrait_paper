#Process genome_scan/summstats
#and collate results into sqlite3 
#database

import glob
import sqlite3
import pandas as pd
import sys,os

opt = sys.argv[1]
mu = sys.argv[2]

#H2_1.0_OPT_0.5_mu_0.001_sigmu0.25_10regions_popgen.db
outfile = '../../sqlite/H2_1.0_OPT_' + str(opt) + '_mu_' + str(mu) + '_sigmu0.25_10regions_popgen.db'

if os.path.exists(outfile):
    os.remove(outfile)

con = sqlite3.connect(outfile)

#These loop are a bit inefficient
#b/c we're globbing out just 1 file
#at a time.  I'll live with that b/c
#it sorts things as I need.

#Contents of each file are as follows:
#10 rows / time point
#Each row is a locus
#Sampled every 100 generations pre-optimum
#shift
#Sampled every 10 generations post-optimum
#shift for 4N generations.

loci=[i for i in range(10)]*int(50000/100)
temp=[[i]*10 for i in range(0,50000,100)]
gen=[item+100 for sublist in temp for item in sublist]
for batch in range(16):
    for rep in range(64):
        pattern_pre = 'H2*OPT_' + str(opt) + '_mu_' + str(mu) + '*.pre_shift.batch' + str(batch) + '.' + str(rep) + '.summstats.gz'
        files = glob.glob(pattern_pre)
        for i in files:
            df = pd.read_csv(i,sep=' ')
            df['locus'] = loci
            df['generation']=gen
            df['opt'] = [float(opt)]*len(df.index)
            df['mu'] = [float(mu)]*len(df.index)
            df['rep'] = [(batch*64) + rep+1]*len(df.index)
            df.to_sql('data',con,if_exists='append',index=False)

loci=[i for i in range(10)]*int(20010/10)
temp=[[i]*10 for i in range(50000,70010,10)]
gen=[item for sublist in temp for item in sublist]
for batch in range(16):
    for rep in range(64):
        pattern_post = 'H2*OPT_' + str(opt) + '_mu_' + str(mu) + '*.post_shift.batch' + str(batch) + '.' + str(rep) + '.summstats.gz'
        files = glob.glob(pattern_post)
        for i in files:
            df = pd.read_csv(i,sep=' ')
            df['locus'] = loci
            df['generation']=gen
            df['opt'] = [float(opt)]*len(df.index)
            df['mu'] = [float(mu)]*len(df.index)
            df['rep'] = [(batch*64) + rep+1]*len(df.index)
            df.to_sql('data',con,if_exists='append',index=False)

con.execute('create index gen on data (generation)')

con.close()
