import pandas as pd
import multiprocessing as mp

LOCK=mp.Lock()

def process_file(args):
    g=args[0]
    fn=args[1]
    opt=args[2]
    mu=args[3]
    OUTFILE=args[4]
    temp=pd.DataFrame(g)
    d=pd.read_hdf(fn)
    d.reset_index(inplace=True)
    d.set_index(['rep','locus'],inplace=True,drop=True)
    temp.set_index(['rep','locus'],inplace=True,drop=True)
    result= temp.join([d],how='inner').reset_index()
    means=result.groupby(['generation','variable']).mean().reset_index()
    means['opt']=[opt]*len(means.index)
    means['mu']=[mu]*len(means.index)
    del means['rep']
    del means['locus']
    LOCK.acquire()
    out=pd.HDFStore(OUTFILE,'a',complevel=6,complib='zlib')
    out.append('fixations',means)
    out.close()
    LOCK.release()
    temp=None
    d=None
    result=None
    means=None

def process_fixations(x,have_fixations = True):
    #This is the expected range of loci.
    #if any are missing in x for a specific rep, 
    #then those loci did not have sweeps from standing varation
    #nor any fixations after ('hard sweeps')
    loci=set(range(10))
    temp=[]
    xg=x.groupby(['opt','mu'])
    for n,g in xg:
        reps=g.groupby(['rep'])
        for n2,g2 in reps:
            #Get loci that DID 
            #have fixations for this rep
            l=set(g2.locus)
            locus_list = None
            if have_fixations is True:
                locus_list = loci.intersection(l)
            else:
                locus_list = loci.difference(l)
            for f in locus_list:
                t = {'opt':n[0],'mu':n[1],'rep':n2,
                    'locus':f}
                temp.append(t)
    rv=pd.DataFrame(temp)
    return rv

def count_soft_hard(x):
    """
    Counts up # sweeps from
    standing varation and from new mutation
    for each replicate w/in each parameter combo
    """
    soft=x[x.type=='s']
    hard=x[x.type=='h']

    agg={'type':'count'}

    sg = soft.groupby(['opt','mu','rep','locus']).agg(agg).reset_index()
    hg = hard.groupby(['opt','mu','rep','locus']).agg(agg).reset_index()

    sg.set_index(['opt','mu','rep','locus'],inplace=True)
    hg.set_index(['opt','mu','rep','locus'],inplace=True)
    sg.rename(columns={'type':'soft'},inplace=True)
    hg.rename(columns={'type':'hard'},inplace=True)
    result=sg.join(hg,how='outer')
    result.reset_index(inplace=True)
    result.fillna(0,inplace=True)
    return result 
