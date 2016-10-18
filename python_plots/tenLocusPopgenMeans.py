import pandas as pd
import glob
from common_functions import *
import multiprocessing as mp
import gc,tarfile,sys

#Step 1: read in the tar file contents,
#collecting file names according to
#input parameter values.

filenames=[]

tfile = tarfile.open("../genome_scan_summstats.tar","r")

for m in tfile.getmembers():
    n=str.split(m.name,"_")
    if(len(n)>2):
        OPT=n[4]
        mu=n[6]
        #Get infer the replicate number from the file name
        temp=str.split(n[10],".")
        batch=int(temp[1].replace("batch",""))
        subrep=int(temp[2])
        rep=batch*64+subrep
        filenames.append( {'opt':OPT,'mu':mu,"fn":m.name,"rep":rep} )

files=pd.DataFrame(filenames)

#Split files up based on before or after optimum shift.
prefiles=files[files.fn.str.contains("pre")]
postfiles=files[files.fn.str.contains("post")]

def processFiles(x,tstart,interval):
    #There are 10 loci.
    #Sampling starts at tstart
    #and continued every "interval"
    #from then on
    df=[]
    I=0
    for xi in g.fn:
        #print xi
        f=tfile.extractfile(xi)
        si=pd.read_csv(f,sep=" ",compression='gzip')
        #print "read",len(si.index),I,len(x.index),x.iloc[I].rep
        ntimepoints = len(si.index)/10
        loci=range(10)*ntimepoints
        generation=[]
        for i in range(ntimepoints):
            generation.extend([tstart+i*interval]*10)
        si['locus']=loci
        si['generation']=generation
        si['rep']=[x.iloc[I].rep]*len(si.index)
        I+=1
        df.append(si)
    return pd.concat(df)
    #print type( df[0])
    #dfc=pd.concat(df)
    #print type(dfc)
    #dfm=pd.melt(dfc,id_vars=['generation','locus','rep'])
    #return dfm.groupby(['generation','variable']).mean().reset_index()

#Group by params
prefilesg = prefiles.groupby(['opt','mu'])
postfilesg = postfiles.groupby(['opt','mu'])

#Now, we go through all the files, 
#and collect a few types of output:
#1. Raw values across time & loci & params, with parameters as index
#2. Mean values across time & params, with raw index
#Data will be stored 'melted'

def updateMeans(x,n,f):
    xm=pd.melt(x,id_vars=['generation','locus','rep'])
    xg=xm.groupby(['generation','variable']).mean().reset_index()
    xg['opt']=[n[0]]*len(xg.index)
    xg['mu']=[n[1]]*len(xg.index)
    xg['sigmu']=[0.25]*len(xg.index)
    xg.drop(['locus'],axis=1,inplace=True)
    f.append('means',xg)

def updateRaw(x,n,f):
    xm=pd.melt(x,id_vars=['generation','locus','rep'])
    xm['opt']=[n[0]]*len(xm.index)
    xm['mu']=[n[1]]*len(xm.index)
    xm['sigmu']=[0.25]*len(xm.index)
    xm.set_index(['opt','mu'],inplace=True,drop=True)
    print xm.head()
    f.append('raw',xm)

meanFile=pd.HDFStore('tenLocusPopgenMeans.h5','w',complevel=6,complib='zlib')
rawFile=pd.HDFStore('tenLocusPopgenRaw.h5','w',complevel=6,complib='zlib')

for n,g in prefilesg:
    print n
    x=processFiles(g,100,100)
    print x.generation.max()
    updateMeans(x,n,meanFile)
    updateRaw(x,n,rawFile)

for n,g in postfilesg:
    x=processFiles(g,50000,10)
    updateMeans(x,n,meanFile)
    updateRaw(x,n,rawFile)

