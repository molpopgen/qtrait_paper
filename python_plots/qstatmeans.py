import pandas as pd
import gc,multiprocessing as mp


#fxn to pull simulation params out of file names
def parse_params(fn):
    x=fn.split('_')
    opt=x[3]
    mu=x[5]
    return {'opt':opt,'mu':mu,'sigmu':0.25}

#def update_totals(filei,Totals,START,STOP):
def update_totals(args):
    filei,Totals,START,STOP=args
    print START,STOP
    #print filei,START,STOP
    x=pd.read_hdf(filei,where='rep>=START&rep<STOP')
    x.reset_index(inplace=True)
    x.drop_duplicates(inplace=True)
    xm=x.groupby(['generation','stat']).mean().reset_index()
    if Totals is None:
        Totals=xm
    else:
        Totals.value+=xm.value
    xm=None
    x=None
    del x
    del xm
    return Totals

def get_values(args):
    filei,START,STOP=args
    print START,STOP
    #print filei,START,STOP
    x=pd.read_hdf(filei,where='rep>=START&rep<STOP')
    x.reset_index(inplace=True)
    x.drop_duplicates(inplace=True)
    xm=x.groupby(['generation','stat']).mean().reset_index()
    n=None
    del x
    return xm.value

LOCK=mp.Lock()
#def process_file(out,filei,pi):
def process_file(args):
    print args
    outfn,filei,pi=args
    START=0
    STOP=64
    Totals=update_totals((filei,None,START,STOP))
    START+=64
    STOP+=64
    while START<1024:
        Totals.value += get_values((filei,START,STOP))
        #Totals=P.map(update_totals,[(filei,Totals,START,STOP)])[0]
        START+=64
        STOP+=64
    for k in pi:
        Totals[k]=[pi[k]]*len(Totals.index)
    Totals.value/=float(1024/64)
    LOCK.acquire()
    out=pd.HDFStore(outfn,'a',complevel=6,complib='zlib')
    out.append('stats',Totals)    
    out.close()
    LOCK.release()
    Totals=None
    print "returning"
