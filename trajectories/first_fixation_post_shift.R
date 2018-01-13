# Get the time to the first fixation after the optimum shift.
# Once we have that, we will get the mean trait value in a
# population up until the mean/median time to the first fixation.

library(readr)
library(dplyr)
library(stringr)

summarise_first_fixation_post_shift <- function(f)
{
    f2 = str_replace(f,"traj_merged_recent_fixations","stats")

    db = DBI::dbConnect(RSQLite::SQLite(),f)
    dbt = tbl(db,'freqs')

    first_post_shift_q = dbt %>% filter(freq == 1,generation > 50000)

    first_post_shift = collect(first_post_shift_q) %>% group_by(repid) %>%
        arrange(generation) %>%
        filter(row_number() == 1) %>%
        select(-index,-freq)

    # Let's refine things and get the first time the first fixation
    # crosses p = 0.1 and p = 0.5

    first_pass_q = dbt %>% filter(generation > 50000,freq >= 0.1)

    first_pass = collect(first_pass_q) %>%
        group_by(repid,origin,pos,esize) %>%
        filter(row_number() == 1)

    j = first_post_shift %>% left_join(first_pass,by=c("repid","pos","esize","origin")) 

    mean10 = mean(j$generation.y)
    median10 = median(j$generation.y)

    first_pass_q = dbt %>% filter(generation > 50000,freq >= 0.5)

    first_pass = collect(first_pass_q) %>%
        group_by(repid,origin,pos,esize) %>%
        filter(row_number() == 1)

    j = first_post_shift %>% left_join(first_pass,by=c("repid","pos","esize","origin")) 

    mean50 = mean(j$generation.y)
    median50 = median(j$generation.y)

    DBI::dbDisconnect(db)

    mean_first_ftime = mean(first_post_shift$generation)
    median_first_ftime = median(first_post_shift$generation)
    trait_db = paste("../sqlite/",f2,sep="")
    db = DBI::dbConnect(RSQLite::SQLite(),trait_db)
    dbt = tbl(db,'data')

    mt = ceiling(max(mean_first_ftime,median_first_ftime))
    trait_value_query = dbt %>% filter(stat=='tbar',
                                       generation >= 50000,
                                       generation <= mt,
                                       rep < 256)

    s = str_split(f,"_")
    opt=as.numeric(s[[1]][4])
    mu=as.numeric(s[[1]][6])
    trait_value_data = collect(trait_value_query) %>% 
        group_by(generation) %>%
        summarise(mz = mean(value)) %>% 
        mutate(mean10=mean10,median10=median10,
                mean50=mean50,median50=median50,
                mean_first_ftime=mean_first_ftime,
                median_first_ftime=median_first_ftime,
                opt=opt,mu=mu)
    DBI::dbDisconnect(db)
    trait_value_data
}

files = dir('.','merged_recent_fixations.db$')
I=0
for (i in files)
{
    x = summarise_first_fixation_post_shift(i)
    APPEND = ifelse(I==0,FALSE,TRUE)
    write_delim(x,"first_fixation_properties.txt",delim=" ",append=APPEND)
    I=I+1
}
