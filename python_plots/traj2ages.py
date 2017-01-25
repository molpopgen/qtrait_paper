import sqlite3
import pandas as pd
import multiprocessing as mp
import sys

def run_query(args):
    dbfile,repid=args
    print repid
    conn=sqlite3.connect(dbfile)
    sql='select generation,count(generation) as nobs,avg(generation-origin) as mean_age from freqs where rep == '+str(repid)+' group by generation;'
    rv= pd.read_sql_query(sql,conn)
    rv['repid']=[repid]*len(rv.index)
    conn.close()
    return rv

if __name__=="__main__":
    dbfile=sys.argv[1]

    P=mp.Pool()
    args=[(dbfile,i) for i in range(200)]
    results=P.map(run_query,args,25)
    P.close()
    P.join()
    allres=[res for res in results]
    print len(allres)
