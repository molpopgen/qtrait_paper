#HDF5 -> .db conversion
#Output data are indexed by (rep,generation),
#and the resulting db has a single table called "data"
import pandas as pd
import sys,sqlite3

if __name__=="__main__":
    infile=sys.argv[1]
    dbname=sys.argv[2]

    x=pd.read_hdf(infile,chunksize=50000,iterator=True)
    conn=sqlite3.connect(dbname)
    for chunk in x:
        chunk.to_sql('data',conn,if_exists='append')
    conn.close()

