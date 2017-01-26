import sqlite3
import pandas as pd
import numpy as np
import multiprocessing as mp
import sys

def run_query(args):
    dbfile,repid=args
    conn=sqlite3.connect(dbfile)
    sql='select generation,count(generation) as nobs,avg(generation-origin) as mean_age from freqs where rep == '+str(repid)+' group by generation;'
    rv= pd.read_sql_query(sql,conn)
    rv['repid']=[repid]*len(rv.index)
    conn.close()
    return rv

if __name__=="__main__":
    dbfile=sys.argv[1]
    outfile=sys.argv[2]
    P=mp.Pool()
    args=[(dbfile,i) for i in range(1024)]
    results=P.map(run_query,args,16)
    P.close()
    P.join()
    allres=pd.concat([res for res in results])
    allres=allres[allres.generation >= 8*5000]
    ttl_obs = allres.groupby(['generation']).sum().reset_index()
    allres.set_index(['generation'],inplace=True,drop=True)
    ttl_obs.set_index(['generation'],inplace=True,drop=True)
    ttl_obs.drop(['repid','mean_age'],axis=1,inplace=True)
    combined=allres.join([ttl_obs],how='left')
    combined['pnobs']=combined.nobs_x/combined.nobs_y
    combined['weighted_mean_age']=combined.mean_age*combined.pnobs
    combined.drop(['repid','nobs_x','nobs_y','pnobs'],axis=1,inplace=True)
    combined.reset_index(inplace=True)
    combined=combined.groupby(['generation']).sum().reset_index()
    combined.drop(['mean_age'],axis=1,inplace=True)
    combined.to_csv(outfile,sep=' ',header=True,index=False,compression='gzip')
