#!/usr/bin/env/Rscript

suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(stringr))
library(RSQLite)
library(methods)
args <- commandArgs(trailing=TRUE)

#make file names for outputs
fixations_db_name = "output/fixations.db"
fixations_db = src_sqlite(fixations_db_name,create=TRUE)
ages_esizes_db_name = "output/ages_esizes.db"
ages_esizes_db = src_sqlite(ages_esizes_db_name,create=TRUE)

collect_fixations <- function(xx,mu,opt)
{
    fixations <- r %>% filter(freq == 1.0) %>% arrange(origin) %>%
        mutate(mu=as.numeric(mu)) %>% mutate(opt=as.numeric(opt))
    fixations
}

collect_ages_esizes <- function(xx,mu,opt)
{
    ages <- xx %>% group_by(generation,repid) %>%
        summarise(mean_age = mean(generation-origin),mesize=mean(esize),nmuts=n()) %>%
        select(generation,mean_age,mesize,nmuts,repid) %>%
        arrange(generation) %>%
        mutate(mu=as.numeric(mu)) %>% mutate(opt=as.numeric(opt))
    ages
}

for (infile in args)
{
    s = str_split(infile,"_")
    opt = s[[1]][4]
    mu=s[[1]][6]
    
    db <- src_sqlite(infile)
    db_tbl <- tbl(db,'freqs')

    #Do things per replicate
    #so that total RAM use is kept sane.
    #(This is slower, though)
    for(id in seq(0,1023))
    {
        q <- db_tbl %>% 
    	    filter(repid == id)
    
        r = collect(q,n=Inf)
    
        #pull out fixations and append them to our database
        fixations <- collect_fixations(r,mu,opt)
        dbWriteTable(fixations_db$con,"fixations",as.data.frame(fixations),append=TRUE)

        #Now, we get allele ages, per generation and per replicate
        all_ages <- collect_ages_esizes(r,mu,opt)

        #Now, we get allele ages and effect sizes,
        #excluding mutations that eventually fix

        #Pull out mutations that did not get fixed,
        notfixed <- r %>% group_by(pos,esize,origin) %>%
            filter(max(freq) < 1.0)

        #Get a table of age/esize data with and w/o fixations
        #via a join
        ages_table = collect_ages_esizes(notfixed,mu,opt) %>% 
            select(generation,mean_age,mesize) %>%
            rename(mean_age_no_fixations = mean_age) %>%
                rename(mesize_no_fixations = mesize) %>%
                left_join(all_ages,by="generation")

        dbWriteTable(ages_esizes_db$con,"ages_esizes",as.data.frame(ages_table),append=TRUE)
    }
}
