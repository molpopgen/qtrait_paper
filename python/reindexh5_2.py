#This file is for the fixations files only
import pandas as pd
import sys

infile = sys.argv[1]
outfile = sys.argv[2]
dfname = sys.argv[3]
fitr = pd.read_hdf(infile,iterator=True,chunksize=500000)
out = pd.HDFStore(outfile,'w',complevel=6,complib='zlib')
for c in fitr:
    c.set_index(['rep'],drop=True,inplace=True)
    out.append(dfname,c)

