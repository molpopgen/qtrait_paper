#!/usr/bin/env r

suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(RSQLite))
suppressMessages(library(docopt))

source("filename2params.R")

doc <- "Usage: summarise_10locus_fixations.R [-i <infile> -o <outfile>]"

options <- docopt(doc)

db <- src_sqlite(options$infile)

data <- tbl(db,'fixations')

#Get data on loci with any sweeps at all
#with respect to optimum shift

non_neutral_fixations = data %>%
    filter(ftime > 50000 & neutral == 0) %>%
    mutate(locus = as.integer(pos))

sweeps = collect(non_neutral_fixations)

hard_sweeps = sweeps %>%
    filter(g < 50000)

soft_sweeps = sweeps %>%
    filter(g > 50000)

#Count no. of each type of sweep per locus 
#and per replicate

hard_sweeps_per_locus = hard_sweeps %>%
    select(locus,rep) %>%
    group_by(locus,rep) %>%
    mutate(nhard=n())

soft_sweeps_per_locus = soft_sweeps %>%
    select(locus,rep) %>%
    group_by(locus,rep) %>%
    mutate(nsoft=n())

#Combine the last two results to get total of each
#type of sweep per locus

ttl_sweeps_per_locus = inner_join(hard_sweeps_per_locus,soft_sweeps_per_locus,by=c('rep','locus')) %>% distinct(rep,locus)

#Still not done: need to fill in zero for loci
#that had no sweeps at all

rep_locus = data.frame(locus=rep(0:9,1024),rep=rep(0:1023,each=10))

all_data = full_join(rep_locus,ttl_sweeps_per_locus,by=c('rep','locus'))

#Joins put NA where there is no overlap,
#which means no sweeps, so we set to 0
all_data$nhard[is.na(all_data$nhard)]=0
all_data$nsoft[is.na(all_data$nsoft)]=0

gzo=gzfile(options$outfile,"w")
write.table(all_data,gzo,quote=F,row.names=F)
close(gzo)
