library(dplyr)
library(readr)

db=DBI::dbConnect(RSQLite::SQLite(),"all_fixation_summaries.sqlite3")
dbt=tbl(db,"fixations")

soft=collect(dbt) %>%
    group_by(opt,mu,repid)%>%
    filter(origin+sojourn_time-1>50000)%>%
    mutate(n=n())%>%
    select(-pos,-origin) %>%
    group_by(opt,mu) %>%
    summarise(mesize_soft=mean(esize),msojourn_soft=mean(sojourn_time)/5e3,mean_soft=mean(n))

hard=collect(dbt) %>%
    group_by(opt,mu,repid)%>%
    filter(origin>50000)%>%
    mutate(n=n())%>%
    select(-pos,-origin) %>%
    group_by(opt,mu) %>%
    summarise(mesize_hard=mean(esize),msojourn_hard=mean(sojourn_time)/5e3,mean_hard=mean(n))

old=collect(dbt) %>%
    group_by(opt,mu,repid)%>%
    filter(origin+sojourn_time-1<50000)%>%
    mutate(n=n())%>%
    select(-pos,-origin) %>%
    group_by(opt,mu) %>%
    summarise(mesize_old=mean(esize),msojourn_old=mean(sojourn_time)/5e3,mean_old=mean(n))

final<-soft
final$mean_hard<-hard$mean_hard
final$msojourn_hard<-hard$msojourn_hard
final$mesize_hard<-hard$mesize_hard

final$mean_old<-old$mean_old
final$msojourn_old<-old$msojourn_old
final$mesize_old<-old$mesize_old
print(final)
write_delim(round(final,digits=4),"traj_fixation_summary_table.txt",delim=" ")
