library(feather)
library(dplyr)

x=read_feather("standing_vs_new.feather")

means = x %>% group_by(opt,rep,type) %>% mutate(nt=length(type)) %>% group_by(opt,type) %>% summarise(mnt = mean(nt),mft=mean(flen),mesize=mean(esize),mfo=mean(f0))

print(means)

medians = x %>% group_by(opt,rep,type) %>% mutate(nt=length(type)) %>% group_by(opt,type) %>% summarise(mnt = median(nt),mft=median(flen),mesize=median(esize),mfo=median(f0))

print(medians)
