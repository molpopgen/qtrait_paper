library(dbplyr)
library(dplyr)
library(tidyr)
library(lattice)
library(viridis)

dbf = src_sqlite("genome_scan_with_sweep_counts.sqlite3")
dbt = tbl(dbf,'data')

hard_rec_only_q = dbt %>%
    filter(h100 >0 & ttl_sweeps == h100) %>%
    group_by(opt,mu,generation,window) %>%
    summarise(hmeanD = mean(tajd), hmeanH = mean(Hprime))

hard_rec_only = collect(hard_rec_only_q,n=Inf) %>%
    gather(key = statistics, value = value, hmeanD, hmeanH)

soft_q = dbt %>%
    filter(soft >0 & ttl_sweeps == soft) %>%
    group_by(opt,mu,generation,window) %>%
    summarise(smeanD = mean(tajd), smeanH = mean(Hprime))

soft = collect(soft_q,n=Inf) %>%
    gather(key = statistics, value = value, smeanD, smeanH)

nosweeps_q = dbt %>%
    filter(ttl_sweeps == 0) %>%
    group_by(opt,mu,generation,window) %>%
    summarise(nosweep_meanD = mean(tajd), nosweep_meanH = mean(Hprime))

nosweeps = collect(nosweeps_q,n=Inf) %>%
    gather(key = statistics, value = value, nosweep_meanD, nosweep_meanH)

# Filtering on window 5 here
# means we focus only where
# selected mutations are happening
data = rbind(hard_rec_only,soft,nosweeps) %>%
    filter(window == 5) %>%
    mutate(scaled_time = (generation - 5e4)/5e4)


KEY=list(space="top",columns=3,title="",
         cex.title=1,points=FALSE,lines=TRUE,just=0.5)
STRIP=strip.custom(strip.names = TRUE, 
                   var.name = c(expression(z[o]),expression(mu)),bg=c("white"))
COLORS=viridis(length(unique(as.factor(data$statistics))))
p = xyplot(value ~ scaled_time|as.factor(opt)*as.factor(mu),data=data,
           type='l',
           lwd=3,
           group = statistics, par.settings=simpleTheme(col=COLORS),
          auto.key=KEY, xlab="Time since optimum shift (units of N generations)",
          ylab="Mean value of statistic.",
          scales=list(cex=1,alternating=F),
          strip=STRIP,xlim=c(-0.05,0.2))

trellis.device(device="png",file="test.png")
print(p)
dev.off()
