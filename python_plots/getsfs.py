#Get SFS for neutral mutations.
#The input is a big .tar file,
#and the ouptut will be a new, big, 
#.tar file.
import gzip
import sys
import tarfile,string,os
from subprocess import call

#This is a .tar file of the .gz files
#that are the time series "ms" output
tarfilename = sys.argv[1]
#this is the tar file where the SFS will go
otarfilename = sys.argv[2]

tf = tarfile.open(tarfilename,'r')
otf = tarfile.open(otarfilename,'w')

#Get each element in the input tar file
tfiles=[]
for f in tf.getmembers():
    #Check that it matches the pattern
    if f.name.find('neutral')>0:
        if f.name.find('gz')>0:
            tfiles.append(f)

#We will process BATCHSIZE files at a time
#using an external C++ program that will
#generate the SFS for each replicate in 
#each file.  The output file name willl be
#auto-created: foo.gz -> foo.sfs.gz.
#Each output file will be appended to a .tar file
#After each "batch", we clean up.
BATCHSIZE=300
for i in range(0,len(tfiles),BATCHSIZE):
    for j in tfiles[i:i+BATCHSIZE]:
        tf.extract(j)
    #create the command line, 'c',
    #which is a call to an external program
    #based on Intel's TBB library to do the
    #hard work that would be slow in Python
    #30 = no. cores to use.
    c=['../src/sfs','30']
    c.extend([j.name for j in tfiles[i:i+BATCHSIZE]])
    call(c)

    #Cleanup
    for j in tfiles[i:i+BATCHSIZE]:
        os.remove(j.name)
        of=str.replace(j.name,"gz","sfs.gz")
        otf.add(of)
        os.remove(of)

tf.close()
otf.close()
