#!/usr/bin/env r

suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(RSQLite))
suppressMessages(library(docopt))
source("filename2params.R")

doc <- "Usage: summarise_qstats.R [-i <infile> -o <outfile>]"

options <- docopt(doc)

db <- src_sqlite(options$infile)

data <- tbl(db,'data')

query <- data %>%
    filter(generation >= 8*5000) %>%
    group_by(generation,stat) %>%
    summarise(mvalue = mean(value))

results = collect(query)

params = getparams(options$infile)
results.tidy = spread(results,key=stat,value=mvalue) %>%
	mutate(opt=params$opt) %>%
	mutate(mu=params$mu)

gzo=gzfile(options$outfile,"w")

write.table(results.tidy,gzo,quote=FALSE,
           row.names=FALSE)

close(gzo)

