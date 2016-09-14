import pandas as pd
import gc,multiprocessing

#fxn to pull simulation params out of file names
def parse_params(fn):
    x=fn.split('_')
    opt=x[3]
    mu=x[5]
    return {'opt':opt,'mu':mu,'sigmu':0.25}

#def update_totals(filei,Totals,START,STOP):
def update_totals(args):
    filei=args[0]
    Totals=args[1]
    START=args[2]
    STOP=args[3]
    print filei,START,STOP
    x=pd.read_hdf(filei,where='rep>=START&rep<STOP')
    x.reset_index(inplace=True)
    x.drop_duplicates(inplace=True)
    xm=x.groupby(['generation','stat']).mean().reset_index()
    if Totals is None:
        Totals=xm
    else:
        Totals.value+=xm.value
    return Totals

def process_file(out,filei,pi):
    START=0
    STOP=64
    Totals=None
    gc.collect()
    while START<1024:
        P=multiprocessing.Pool(1)
        Totals=P.map(update_totals,[(filei,Totals,START,STOP)])[0]
        P.close()
        START+=64
        STOP+=64
    for k in pi:
        Totals[k]=[pi[k]]*len(Totals.index)
    Totals.value/=float(1024/64)
    out.append('stats',Totals)    
