library(feather)
library(dplyr)

x=read_feather("standing_vs_new.feather")

meanTypes = x %>% group_by(opt,rep,type) %>% mutate(nt=length(type)) %>% group_by(rep,nt,opt,type) %>% summarise(foo=length(nt)) %>% group_by(opt,type) %>% summarise(mtype=mean(foo),mintype=min(foo),maxtype=max(foo))

print(meanTypes)
