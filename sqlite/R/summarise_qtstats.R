#!/usr/bin/env r

suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(RSQLite))
suppressMessages(library(docopt))
source("filename2params.R")

doc <- "Usage: summarise_loads.R [-i <infile> -o <outfile>]"

opt <- docopt(doc)

db <- src_sqlite(opt$infile)

data <- tbl(db,'data')

query <- data %>%
    filter(generation >= 8*5000) %>%
    group_by(generation,stat) %>%
    summarise(mvalue = mean(value))

results = collect(query,n=5)

params = getparams(opt$infile)
results.tidy = spread(results,key=stat,value=mvalue) %>%
	mutate(opt=params$opt) %>%
	mutate(mu=params$mu)

gzo=gzfile(opt$outfile,"w")

write.table(results.tidy,gzo,quote=FALSE,
           row.names=FALSE)

close(gzo)

