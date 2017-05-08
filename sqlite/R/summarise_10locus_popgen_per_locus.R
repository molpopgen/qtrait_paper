#!/usr/bin/env r

suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(RSQLite))
suppressMessages(library(docopt))
source("filename2params.R")

doc <- "Usage: summarise_10locus_popgen_per_locus.R [-i <infile> -o <outfile>]"

options <- docopt(doc)

db <- src_sqlite(options$infile)

data <- tbl(db,'data')

query <- data %>%
    filter(generation >= 8*5000) %>%
    group_by(generation,locus,variable) %>%
    summarise(mvalue = mean(value))

results = collect(query)

params = getparams(options$infile)
results.tidy = spread(results,key=variable,value=mvalue) %>%
	mutate(opt=params$opt) %>%
	mutate(mu=params$mu)

gzo=gzfile(options$outfile,"w")

write.table(results.tidy,gzo,quote=FALSE,
           row.names=FALSE)

close(gzo)

