library(knitr)
library(dplyr)
library(readr)
library(ggplot2)

db=DBI::dbConnect(RSQLite::SQLite(),"all_fixation_summaries.sqlite3")
dbt=tbl(db,"fixations")

x = collect(dbt) %>%
    mutate(scaled = 5e3*esize*esize,
           Category = ifelse(origin < 5e4,"Standing variation",
                             ifelse(origin<=5e4+100,"New mutation","None"))) %>%
    filter(Category != "None", scaled <= 1e3) %>%
    group_by(opt,mu,repid,strength=cut(scaled,c(0, 10, 100, 1000, 1e8)),Category) %>%
    mutate(n=n()) %>%
    group_by(opt, mu, strength, Category) %>%
    summarise(mean_sweeps = mean(n),msojourn = mean(sojourn_time)/5e3)

xlabels=c("(0,10]","(10,100]","(100,1000]",">= 1000")

number = ggplot(x, aes(x=strength,y=mean_sweeps)) +
    facet_grid(opt ~ mu,
               labeller=label_bquote(rows=z[o]~"="~.(opt),cols=mu~"="~.(mu))) +
    geom_point(alpha=0.5,aes(color=Category)) +
    xlab(expression(paste('N',gamma^2,' range'))) +
    ylab("Mean number of fixations") +
    scale_x_discrete(labels=xlabels) + theme_bw() +
    theme(axis.text.x = element_text(angle=45,hjust=1)) #+
ggsave('NumFixationsSingleRegion.pdf')


number = ggplot(x, aes(x=strength,y=msojourn)) +
    facet_grid(opt ~ mu,
               labeller=label_bquote(rows=z[o]~"="~.(opt),cols=mu~"="~.(mu))) +
    geom_point(alpha=0.5,aes(color=Category)) +
    xlab(expression(paste('N',gamma^2,' range'))) +
    ylab("Mean fixation time (units of N generations)") +
    scale_x_discrete(labels=xlabels) + theme_bw() +
    theme(axis.text.x = element_text(angle=45,hjust=1)) #+
ggsave('FixationTimesSingleRegion.pdf')
q('no')

soft=collect(dbt) %>%
    mutate(large = ifelse(abs(esize) >= 2*sqrt(2)*sqrt(mu),T,F)) %>%
    group_by(opt,mu,repid,large)%>%
    filter(origin+sojourn_time-1>50000)%>%
    mutate(n=n())%>%
    select(-pos,-origin) %>%
    group_by(mu,opt,large) %>%
    summarise(mesize_soft=mean(abs(esize)),msojourn_soft=mean(sojourn_time)/5e3,mean_soft=mean(n))

hard=collect(dbt) %>%
    mutate(large = ifelse(abs(esize) >= 2*sqrt(2)*sqrt(mu),T,F)) %>%
    group_by(opt,mu,repid,large)%>%
    filter(origin>50000 & origin < 5e4 + 1e2)%>%
    mutate(n=n())%>%
    select(-pos,-origin) %>%
    group_by(mu,opt,large) %>%
    summarise(mesize_hard=mean(abs(esize)),msojourn_hard=mean(sojourn_time)/5e3,mean_hard=mean(n))

old=collect(dbt) %>%
    mutate(large = ifelse(abs(esize) >= 2*sqrt(2)*sqrt(mu),T,F)) %>%
    group_by(opt,mu,repid,large)%>%
    filter(origin+sojourn_time-1<50000)%>%
    mutate(n=n())%>%
    select(-pos,-origin) %>%
    group_by(opt,mu,large) %>%
    summarise(mesize_old=mean(abs(esize)),msojourn_old=mean(sojourn_time)/5e3,mean_old=mean(n))

final = left_join(soft,hard,by=c("opt","mu","large")) %>%
    arrange(mu,opt,large)
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
