#HDF5 -> .db conversion
#Output data are indexed by (rep,generation),
#and the resulting db has a single table called "data"
import pandas as pd
import sys,sqlite3

if __name__=="__main__":
    infile=sys.argv[1]
    dbname=sys.argv[2]

    x=pd.read_hdf(infile,chunksize=1000000,iterator=True)
    conn=sqlite3.connect(dbname)
    for chunk in x:
        chunk.reset_index(inplace=True,drop=True)
        indexes=[]
        if 'rep' in chunk:
            indexes.append('rep')
        if 'generation' in chunk:
            indexes.append('generation')
        haveIndex=False
        if len(indexes) > 0:
            chunk.set_index(indexes,inplace=True,drop=True)
            haveIndex=True
        chunk.rename(columns={"total_Aa":"total_het","total_aa":"total_hom"},inplace=True)
        chunk.to_sql('data',conn,if_exists='append',index=haveIndex)
    conn.close()

