# Division by 10 goes from mean
# per replicate to mean per locus
library(readr)
library(dplyr)
library(lattice)
library(ggplot2)

x = read_delim("raw_fixations.txt.gz",delim=" ")

N = 5e3


x = x %>%
    mutate(large = ifelse(abs(s) > 2*sqrt(2)*sqrt(mu),as.integer(1),as.integer(0)))

num_soft_sweeps = x %>%
    filter(sweep_type == 'soft') %>%
    group_by(mu, opt, repid, large) %>%
    mutate(nsoft = n())

mean_soft_sweeps_per_locus = num_soft_sweeps %>%
    group_by(mu,opt,large) %>%
    summarise(mean_soft = mean(nsoft)/10, mfixtime_soft = median(sojourn_time)/N)

num_hard_sweeps = x %>%
    filter(sweep_type == 'hard' & g <= 50100) %>%
    group_by(mu, opt, repid, large) %>%
    mutate(nhard = n())

mean_hard_sweeps_per_locus = num_hard_sweeps %>%
    group_by(mu,opt,large) %>%
    summarise(mean_hard_recent = mean(nhard)/10, mfixtime_recent = median(sojourn_time)/N)

num_hard_sweeps_late = x %>%
    filter(sweep_type == 'hard' & g > 50100) %>%
    group_by(mu, opt, repid, large) %>%
    mutate(nhard = n())

mean_hard_sweeps_per_locus_late = num_hard_sweeps_late %>%
    group_by(mu,opt,large) %>%
    summarise(mean_hard_late = mean(nhard)/10, mfixtime_late = median(sojourn_time)/N)

joined = mean_soft_sweeps_per_locus %>% left_join(mean_hard_sweeps_per_locus,by=c("mu","opt","large")) %>%
    left_join(mean_hard_sweeps_per_locus_late,by=c("mu","opt","large"))

joined[is.na(joined)] <- 0

write_delim(round(joined,digits=3),"mean_numbers_of_sweeps.txt",delim=" ")

    
