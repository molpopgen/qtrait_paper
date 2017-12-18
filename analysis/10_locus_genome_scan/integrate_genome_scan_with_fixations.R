# Notes: sweep_type = none
# means fixation prior to optimum shift

# This script takes the pop gen statistic data 
# over time and tacks on how many sweeps of various 
# types happened at each locus.  The sweeps are 
# categorized as follows: 
# h100 = hard, within first 100 gens of optimum shift
# h200 = hard, > 100 and <= 200 gens post-shift
# h = hard, > 200 gens post-shift
# soft = from standing variation
# ttl_sweeps = sum of all of the above
library(readr)
library(dplyr)
library(tidyr)
library(dbplyr)
library(stringr)
library(viridis)
library(lattice)
library(RSQLite)

process_genome_scan <- function(filename, fixations)
{
    fs <- str_split(filename,regex('\\.'))[[1]]
    fs[2]<-str_replace(fs[2],'mu','')
    fs[4]<-str_replace(fs[4],'opt','')
    dbmu <- as.numeric(paste(fs[2],fs[3],sep='.'))
    dbopt <- as.numeric(paste(fs[4],fs[5],sep='.'))

    f = fixations %>% filter(opt==dbopt,mu==dbmu)
    dbname = paste("../../mlocus_pickle/",filename,sep="")
    db <- src_sqlite(dbname)
    dbt <- tbl(db,'data')

    dt = collect(dbt) %>%
        full_join(f,by=c('repid','locus')) %>%
        # The full join will put NA values
        # for all loci that had no sweeps at all.
        # We need to recode those with 0.
        # The join also mangles mu/opt for loci w/no sweeps,
        # so we need to put that back in.
        replace_na(list(h100=0,h200=0,h=0,soft=0,ttl_sweeps=0,opt = dbopt, mu = dbmu)) %>%
        # replace_na changes our int to numeric
        mutate(h100 = as.integer(h100),h200 = as.integer(h200), h = as.integer(h),
               soft = as.integer(soft), ttl_sweeps = as.integer(ttl_sweeps))
    dt
}

fixations = read_delim("raw_fixations.txt.gz",delim=" ")

#See note in file header above re: why we remove none here
fixations = fixations %>% filter(sweep_type != 'none') %>%
    group_by(repid,locus,mu,opt) %>%
    summarise(h100 = as.integer(length(which(sweep_type == 'hard' & g <= 5e4 + 100))),
              h200 = as.integer(length(which(sweep_type == 'hard' & g > 5e4 + 100 & g <= 5e4 + 200))),
              h = as.integer(length(which(sweep_type == 'hard' & g > 5e4 + 200))),
              soft = as.integer(length(which(sweep_type == 'soft')))) %>%
    mutate(ttl_sweeps = h100 + h200 + h + soft)

files <- dir(path="../../mlocus_pickle/",pattern="*.genome_scan.db")


DBname="genome_scan_with_sweep_counts.sqlite3"
if(file.exists(DBname))
{
    file.remove(DBname)
}
DB = src_sqlite(DBname,create=T)
for (i in files)
{
    dfi = process_genome_scan(i,fixations)
    dbWriteTable(DB$con, "data", dfi,append=TRUE)
}

# create DB indexes
dbSendQuery(DB$con, "create INDEX h100i on data(h100)")
dbSendQuery(DB$con, "create INDEX h200i on data(h200)")
dbSendQuery(DB$con, "create INDEX hi on data(h)")
dbSendQuery(DB$con, "create INDEX si on data(soft)")
