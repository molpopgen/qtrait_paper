library(dplyr)
library(tibble)
library(readr)
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

n = commandArgs(trailing=TRUE)

infile = n[1]
outfile = n[2]

indb = src_sqlite(infile)
indb_table = tbl(indb,'data')

raw_data = collect(indb_table)
power_query = raw_data %>% group_by(generation,locus,window) %>%
    summarise(tajd05 = sum(tajd <= critical_vals$tajd_lower05)/n(),
           tajd01 = sum(tajd <= critical_vals$tajd_lower01)/n(),
           tajd001 = sum(tajd <= critical_vals$tajd_lower001)/n(),
           hprime05 = sum(hprime <= critical_vals$hprime_lower05)/n(),
           hprime01 = sum(hprime <= critical_vals$hprime_lower01)/n(),
           hprime001 = sum(hprime <= critical_vals$hprime_lower001)/n(),
           max_abs_nSL05 = sum(max_abx_nSL <= critical_vals$max_abs_nSL_lower05)/n(),
           max_abs_nSL01 = sum(max_abx_nSL <= critical_vals$max_abs_nSL_lower01)/n(),
           max_abs_nSL001 = sum(max_abx_nSL <= critical_vals$max_abs_nSL_lower001)/n())


power_results = collect(power_query)
write_delim(power_results,outfile,delim="\t")
