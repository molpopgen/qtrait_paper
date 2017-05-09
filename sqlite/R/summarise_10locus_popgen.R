#!/usr/bin/env r

suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(RSQLite))
suppressMessages(library(docopt))
source("filename2params.R")

doc <- "Usage: summarise_10locus_popgen.R [-i <infile> -o <outfile>]"

options <- docopt(doc)

db <- src_sqlite(options$infile)

data <- tbl(db,'data')

params = getparams(options$infile)
query <- data %>%
    #Drop the locus column, as we don't need it
    select(-locus,-rep) %>% 
    filter(generation >= 8*5000) %>%
	group_by(generation) %>%
	summarise_each(funs(mean),-generation,-opt,-mu) %>%
	mutate(opt = params$opt) %>%
	mutate(mu = params$mu)

results = collect(query)
gzo=gzfile(options$outfile,"w")
write.table(results,gzo,quote=FALSE,
           row.names=FALSE)
close(gzo)
