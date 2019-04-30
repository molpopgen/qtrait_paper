# Division by 10 goes from mean
# per replicate to mean per locus
library(readr)
library(dplyr)
library(lattice)
library(ggplot2)

x = read_delim("raw_fixations.txt.gz",delim=" ")

N = 5e3


x = x %>%
    mutate(scaled = 5e3*s*s) %>% filter(scaled <= 1000)

soft_weak = x %>% filter(sweep_type == 'soft', scaled > 0, scaled < 10) %>%
    group_by(mu, opt, repid) %>%
    mutate(nsoft_weak = n())

mean_soft_weak_per_locus = soft_weak %>%
    group_by(mu, opt) %>%
    summarise(mean_soft_weak = mean(nsoft_weak)/10, mfixtime_soft_weak=mean(sojourn_time)/N) %>%
    mutate(upper=10)

y = x %>% filter(sweep_type == 'soft'|(sweep_type=='hard'&g<5e4+100)) %>%
    group_by(sweep_type, mu, opt, repid, r= cut(scaled,c(0, 10, 100, 1e3, 1e8))) %>%
    mutate(nsweeps = n())

z = y %>% group_by(mu, opt, sweep_type,r) %>%
    summarise(msweeps = mean(nsweeps)/10, mtime = mean(sojourn_time)/N)
print(z)
write_delim(z[with(z,order(mu, opt, sweep_type)),],"mean_numbers_of_sweeps.txt",delim=" ", col_names=F)
xlabels=c("(0,10]","(10,100]","(100,1000]",">= 1000")
number = ggplot(z, aes(x=r,y=msweeps)) +
    facet_grid(opt ~ mu,
               labeller=label_bquote(rows=z[o]~"="~.(opt),cols=mu~"="~.(mu))) +
    geom_point(alpha=0.5,aes(color=sweep_type)) +
    xlab(expression(paste('N',gamma^2,' range'))) +
    ylab("Mean number of fixations per locus") +
    scale_x_discrete(labels=xlabels) + theme_bw() +
    theme(axis.text.x = element_text(angle=45,hjust=1)) +
    scale_color_discrete(name="Category",labels=c("New mutation","Standing variation"))
ggsave("MeanNumberOfFixationsPerLocus.pdf")
time = ggplot(z, aes(x=r,y=mtime)) +  
    facet_grid(opt ~ mu,
               labeller=label_bquote(rows=z[w]~"="~.(opt),cols=mu~"="~.(mu))) +
    geom_point(alpha=0.5,aes(color=sweep_type)) +
    xlab(expression(paste('N',gamma^2,' range'))) +
    ylab("Mean fixation time (units of N generations)") +
    scale_x_discrete(labels=xlabels) + theme_bw() +
    theme(axis.text.x = element_text(angle=45,hjust=1)) +
    scale_color_discrete(name="Category",labels=c("New mutation","Standing variation"))
ggsave("MeanFixationTimeMultilocus.pdf")


# 
# num_soft_sweeps = x %>%
#     filter(sweep_type == 'soft') %>%
#     group_by(mu, opt, repid, large) %>%
#     mutate(nsoft = n())
# 
# mean_soft_sweeps_per_locus = num_soft_sweeps %>%
#     group_by(mu,opt,large) %>%
#     summarise(mean_soft = mean(nsoft)/10, mfixtime_soft = median(sojourn_time)/N)
# 
# num_hard_sweeps = x %>%
#     filter(sweep_type == 'hard' & g <= 50100) %>%
#     group_by(mu, opt, repid, large) %>%
#     mutate(nhard = n())
# 
# mean_hard_sweeps_per_locus = num_hard_sweeps %>%
#     group_by(mu,opt,large) %>%
#     summarise(mean_hard_recent = mean(nhard)/10, mfixtime_recent = median(sojourn_time)/N)
# 
# num_hard_sweeps_late = x %>%
#     filter(sweep_type == 'hard' & g > 50100) %>%
#     group_by(mu, opt, repid, large) %>%
#     mutate(nhard = n())
# 
# mean_hard_sweeps_per_locus_late = num_hard_sweeps_late %>%
#     group_by(mu,opt,large) %>%
#     summarise(mean_hard_late = mean(nhard)/10, mfixtime_late = median(sojourn_time)/N)
# 
# joined = mean_soft_sweeps_per_locus %>% left_join(mean_hard_sweeps_per_locus,by=c("mu","opt","large")) %>%
#     left_join(mean_hard_sweeps_per_locus_late,by=c("mu","opt","large"))
# 
# joined[is.na(joined)] <- 0
# 
# write_delim(round(joined,digits=3),"mean_numbers_of_sweeps.txt",delim=" ")
# 
#     
