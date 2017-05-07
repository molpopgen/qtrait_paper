#!/usr/bin/env r

suppressMessages(library(dplyr))
suppressMessages(library(RSQLite))
suppressMessages(library(docopt))

doc <- "Usage: summarise_loads.R [-i <infile> -o <outfile>]"

opt <- docopt(doc)

#if (opt$deps == "TRUE" || opt$deps == "FALSE") {
#    opt$deps <- as.logical(opt$deps)
#} else if (opt$deps == "NA") {
#    opt$deps <- NA
#}

db <- src_sqlite(opt$infile)

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

results = collect(query,n=5)

gzo=gzfile(opt$outfile)

write.table(results,gzo,quote=FALSE,
           row.names=FALSE)

close(gzo)

