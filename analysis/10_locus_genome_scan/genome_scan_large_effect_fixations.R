library(readr)
library(dplyr)
library(tidyr)
library(dbplyr)
library(stringr)
library(viridis)
library(lattice)
library(RSQLite)
library(latticeExtra)

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

    gammahat = 2*sqrt(2)*sqrt(dbmu)
    dt = collect(dbt) %>%
        full_join(f,by=c('repid','locus'))  %>%
        na.omit() 
    dt
}

fixations_raw = read_delim("raw_fixations.txt.gz",delim=" ")

fixations = fixations_raw %>%
    filter(sweep_type != 'none',(g-5e4)/5e4<1e2/5e4,
           abs(s) >= 2*sqrt(2)*sqrt(mu)) %>%
    group_by(repid,locus,mu,opt) %>%
    summarise(nhard=length(which(sweep_type=='hard')),nsoft=length(which(sweep_type=='soft')))


# print(fixations)
# print(subset(fixations,nhard>1&nsoft>1))
# q("no")
# 
# fixations = fixations_raw%>% filter(sweep_type == 'hard') %>% 
#     filter(abs(s) >= 2*sqrt(2)*sqrt(mu)) %>%
#     group_by(repid,locus,mu,opt) %>%
#     summarise(n=n())
# 
# fixations_soft = fixations_raw %>% 
#     filter(sweep_type == 'soft') %>% 
#     filter(abs(s) >= 2*sqrt(2)*sqrt(mu)) %>%
#     group_by(repid,locus,mu,opt) %>%
#     summarise(n=n())
# print(fixations)
# print(fixations_soft)
# q("no")

files <- dir(path="../../mlocus_pickle/",pattern="*.genome_scan.db")


DBnames=c("genome_scan_with_large_effect_hard_sweep_counts.sqlite3",
"genome_scan_with_large_effect_soft_sweep_counts.sqlite3")
DBname="genome_scan_with_large_effect_sweep_counts.sqlite3"
# I=1
# for(DBname in DBnames)
# {
    print(DBname)

    if(file.exists(DBname))
    {
        file.remove(DBname)
    }
    DB = src_sqlite(DBname,create=T)
    for (i in files)
    {
        dfi = process_genome_scan(i,fixations)
        dbWriteTable(DB$con, "data", dfi,append=TRUE)

        # if(I==1)
        # {
        #     dfi = process_genome_scan(i,fixations)
        #     dbWriteTable(DB$con, "data", dfi,append=TRUE)
        # }
        # else
        # {
        #     dfi = process_genome_scan(i,fixations_soft)
        #     dbWriteTable(DB$con, "data", dfi,append=TRUE)
        # }
    }
#    I=I+1
#}
# create DB indexes
# dbSendQuery(DB$con, "create INDEX h100i on data(h100)")
# dbSendQuery(DB$con, "create INDEX h200i on data(h200)")
# dbSendQuery(DB$con, "create INDEX hi on data(h)")
# dbSendQuery(DB$con, "create INDEX si on data(soft)")

