#!/usr/bin/env Rscript

suppressMessages(library(docopt))

doc <- "Usage: combine_single_locus_loads.R [-o <outfile>] [INPUT ...]"	

options <- docopt(doc)

out=gzfile(options$outfile,"wa")
NFILES=0
for (i in options$INPUT)
{
	x=read.table(i,header=T)
    if(NFILES==0)
    {
        write.table(x,out,quote=F,row.names=F,append=FALSE)
    } else {
        write.table(x,out,quote=F,row.names=F,append=TRUE,col.names=FALSE)
    }
    NFILES=NFILES+1
}

close(out)
