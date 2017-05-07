#!/usr/bin/env r

suppressMessages(library(dplyr))
suppressMessages(library(RSQLite))
suppressMessages(library(docopt))
source("filename2params.R")
doc <- "Usage: summarise_loads.R [-i <infile> -o <outfile>]"

options <- docopt(doc)

db <- src_sqlite(options$infile)

data <- tbl(db,'data')

#The columns are perhaps confusing.
#The following are the DB columns, 
#which are means/gen
#total = deviation from fitness of 1
#fixed = load from w calc. only with fixations
#seg = load from w calc. only w/polymorphisms
#total_het = number Aa sites / dip
#total_hom = number aa sites / dip
#total_muts = number mutations/dip

#Our query aggregates the 
#means over sim replicate
query <- data %>%
    filter(generation >= 8*5000) %>%
    group_by(generation) %>%
    summarise(fixed_load = mean(fixed),
              total_load = mean(total),
              seg_load = mean(seg),
              het_sites_per_dip = mean(total_het),
              hom_per_dip = mean(total_hom),
              muts_per_dip = mean(total_muts)
              )
params=getparams(options$infile)
results = collect(query,n=5) %>% 
	mutate(opt = params$opt) %>%
	mutate(mu = params$mu)

gzo=gzfile(options$outfile)

write.table(results,gzo,quote=FALSE,
           row.names=FALSE)

close(gzo)

