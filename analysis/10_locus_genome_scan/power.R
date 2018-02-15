# This script begins the process of calculating
# power per window.
# power.R
# power_gather.R
# plot_unconditional_power.R
# The result is the probability of detecting
# a window at a specific alpha
library(dplyr)
library(tibble)
library(readr)
library(stringr)
ndist_db = src_sqlite("nulldist.db",create=F)
ndist_data_table = tbl(ndist_db, 'data')

lower05 <- function(x)
{
    return(as.numeric(quantile(x,0.05)))
}
lower01 <- function(x)
{
    return(as.numeric(quantile(x,0.01)))
}
lower001 <- function(x)
{
    return(as.numeric(quantile(x,0.001)))
}
raw_dist = collect(ndist_data_table)
critical_vals = raw_dist %>% summarise_all(funs(lower05,lower01,lower001))

files <- dir(path="../../mlocus_pickle/",pattern="*.genome_scan.db")

data = data.frame()
for (infile in files)
{
    print(infile)
    dbname = paste("../../mlocus_pickle/",infile,sep="")
    indb = src_sqlite(dbname)
    indb_table = tbl(indb,'data')

    fs <- str_split(infile,regex('\\.'))[[1]]
    fs[2]<-str_replace(fs[2],'mu','')
    fs[4]<-str_replace(fs[4],'opt','')
    mu <- as.numeric(paste(fs[2],fs[3],sep='.'))
    opt <- as.numeric(paste(fs[4],fs[5],sep='.'))

    raw_data = collect(indb_table)
    power_query = raw_data %>% group_by(generation,locus,window) %>%
        summarise(tajd05 = sum(tajd <= critical_vals$tajd_lower05),
               tajd01 = sum(tajd <= critical_vals$tajd_lower01),
               tajd001 = sum(tajd <= critical_vals$tajd_lower001),
               hprime05 = sum(hprime <= critical_vals$hprime_lower05),
               hprime01 = sum(hprime <= critical_vals$hprime_lower01),
               hprime001 = sum(hprime <= critical_vals$hprime_lower001),
               max_abs_nSL05 = sum(max_abx_nSL <= critical_vals$max_abs_nSL_lower05),
               max_abs_nSL01 = sum(max_abx_nSL <= critical_vals$max_abs_nSL_lower01),
               max_abs_nSL001 = sum(max_abx_nSL <= critical_vals$max_abs_nSL_lower001)) %>%
        mutate(mu=mu,opt=opt)


    power_results = collect(power_query)
    data = bind_rows(data,power_results)
}

of = gzfile("nsig_per_window.txt.gz","w")
write_delim(data,of)
close(of)
