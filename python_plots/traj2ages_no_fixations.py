import sqlite3
import pandas as pd
import numpy as np
import multiprocessing as mp
import sys

def run_query(args):
    dbfile,repid=args
    conn=sqlite3.connect(dbfile)
    sql='select pos,esize,max(freq) as mfreq from freqs where rep == '+str(repid)+' group by pos,esize'
    maxfreqs=pd.read_sql_query(sql,conn)
    sql='select pos,esize,generation,generation-origin as age from freqs where rep == '+str(repid)+' group by pos,esize,generation;'
    ages= pd.read_sql_query(sql,conn)
    maxfreqs.set_index(['pos','esize'],inplace=True)
    ages.set_index(['pos','esize'],inplace=True)
    joined=maxfreqs.join([ages],how='left').reset_index()
    joined=joined[joined.mfreq<1.]
    rv=joined.groupby(['generation']).agg({'age':{'nobs':'count','mean_age':np.mean}}).reset_index()
    rv.columns=rv.columns.droplevel(0)
    rv.columns.values[0]='generation'
    rv['repid']=[repid]*len(rv.index)
    conn.close()
    return rv

if __name__=="__main__":
    dbfile=sys.argv[1]
    outfile=sys.argv[2]
    P=mp.Pool()
    args=[(dbfile,i) for i in range(1024)]
    results=P.map(run_query,args,4)
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
