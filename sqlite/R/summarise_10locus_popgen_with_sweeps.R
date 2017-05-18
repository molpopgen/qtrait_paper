#!/usr/bin/env r

suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(RSQLite))
suppressMessages(library(stringr))
suppressMessages(library(docopt))
source("filename2params.R")

join_and_summarise= function(sweepdata,popgen_data)
{
	rv =semi_join(popgen_data,sweepdata,by=c('rep','locus')) %>%
		select(-locus,-rep) %>%
		group_by(generation) %>%
		summarise_each(funs(mean),-generation)
	rv
}

doc <- "Usage: summarise_10locus_popgen.R [-p <popgen> -s <sweeps>]"

options <- docopt(doc)

#Stub for output
base = options$popgen
base=str_replace(str_replace(base,"../",""),".db","")[[1]]

sweeps = read.table(options$sweeps,header=T)

popgen_db = src_sqlite(options$popgen)

popgen = tbl(popgen_db,'data')

popgen_query = popgen %>% filter(generation >= 40000)

popgen_data = collect(popgen_query)

params = getparams(options$popgen)
nosweeps = join_and_summarise(subset(subset(sweeps,nhard==0),nsoft==0),popgen_data)
hardsweeps = join_and_summarise(subset(subset(sweeps,nhard>0),nsoft==0),popgen_data)
softsweeps = join_and_summarise(subset(subset(sweeps,nhard==0),nsoft>0),popgen_data)

nosweeps_out = paste(base,"_nosweeps.gz",sep="")
gzo=gzfile(nosweeps_out,"w")
write.table(nosweeps,gzo,quote=FALSE,
           row.names=FALSE)
close(gzo)

hardsweeps_out = paste(base,"_hardsweeps.gz",sep="")
gzo=gzfile(hardsweeps_out,"w")
write.table(hardsweeps,gzo,quote=FALSE,
           row.names=FALSE)
close(gzo)

softsweeps_out = paste(base,"_softsweeps.gz",sep="")
gzo=gzfile(softsweeps_out,"w")
write.table(softsweeps,gzo,quote=FALSE,
           row.names=FALSE)
close(gzo)
