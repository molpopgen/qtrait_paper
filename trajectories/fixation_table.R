library(knitr)
library(dplyr)
library(readr)

db=DBI::dbConnect(RSQLite::SQLite(),"all_fixation_summaries.sqlite3")
dbt=tbl(db,"fixations")

soft=collect(dbt) %>%
    mutate(large = ifelse(abs(esize) >= 2*sqrt(2)*sqrt(mu),T,F)) %>%
    group_by(opt,mu,repid,large)%>%
    filter(origin+sojourn_time-1>50000)%>%
    mutate(n=n())%>%
    select(-pos,-origin) %>%
    group_by(opt,mu,large) %>%
    summarise(mesize_soft=mean(abs(esize)),msojourn_soft=mean(sojourn_time)/5e3,mean_soft=mean(n))

hard=collect(dbt) %>%
    mutate(large = ifelse(abs(esize) >= 2*sqrt(2)*sqrt(mu),T,F)) %>%
    group_by(opt,mu,repid,large)%>%
    filter(origin>50000 & origin < 5e4 + 1e2)%>%
    mutate(n=n())%>%
    select(-pos,-origin) %>%
    group_by(opt,mu,large) %>%
    summarise(mesize_hard=mean(abs(esize)),msojourn_hard=mean(sojourn_time)/5e3,mean_hard=mean(n))

old=collect(dbt) %>%
    mutate(large = ifelse(abs(esize) >= 2*sqrt(2)*sqrt(mu),T,F)) %>%
    group_by(opt,mu,repid,large)%>%
    filter(origin+sojourn_time-1<50000)%>%
    mutate(n=n())%>%
    select(-pos,-origin) %>%
    group_by(opt,mu,large) %>%
    summarise(mesize_old=mean(abs(esize)),msojourn_old=mean(sojourn_time)/5e3,mean_old=mean(n))

final = left_join(soft,hard,by=c("opt","mu","large"))
# final = left_join(final,old,by=c("opt","mu","large"))
final$large = as.integer(final$large)
# final<-soft
# final$mean_hard<-hard$mean_hard
# final$msojourn_hard<-hard$msojourn_hard
# final$mesize_hard<-hard$mesize_hard
# 
# final$mean_old<-old$mean_old
# final$msojourn_old<-old$msojourn_old
# final$mesize_old<-old$mesize_old
# print(final)
latex_table = kable(final, format="latex",format.args=list(digits=4))
# write_delim(round(final,digits=4),"traj_fixation_summary_table.txt",delim=" ")
write(latex_table,"traj_fixation_summary_table.txt")
