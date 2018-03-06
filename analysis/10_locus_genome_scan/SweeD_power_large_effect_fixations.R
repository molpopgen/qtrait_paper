library(dplyr)
library(readr)
library(lattice)
library(viridis)

sweep_count_db=src_sqlite('genome_scan_with_large_effect_sweep_counts.sqlite3')
sweep_count_table=tbl(sweep_count_db,'data')

sweep_counts_q = sweep_count_table %>%
    group_by(locus,repid,generation,opt,mu) %>%
    select(nhard,nsoft)

sweep_counts=collect(sweep_counts_q)%>%
    distinct()

sweed=read_delim("SweeD_Sig.txt.gz",delim=" ")
print(unique(sweed$sig05))

sweed_plus_fixations=sweed %>%
    left_join(sweep_counts,by=c('repid','locus','generation','opt','mu')) %>%
    na.omit() #Remove anything w/o large effect fixation

# print(head(sweed_plus_fixations))
hard_only = sweed_plus_fixations %>%
    filter(nhard>0,nsoft==0) %>%
    mutate(sweep_type='New mutation')

soft_only = sweed_plus_fixations %>%
    filter(nhard==0,nsoft>0)%>%
    mutate(sweep_type='Standing variation')

meanLRdata=rbind(hard_only,soft_only) %>%
    group_by(generation,opt,mu,sweep_type) %>%
    summarise(mLR=mean(LR),s05=sum(sig05),n=n()) %>%
    mutate(scaled_time=(generation-5e4)/5e3)

print(meanLRdata[which(meanLRdata$generation == 50000),])
print(unique(meanLRdata$sweep_type))
#COLORS=viridis(length(unique(as.factor(sweed_plus_fixations$sweep_type))))
KEY=list(space="top",columns=2,title="Type of sweep",
         cex.title=1,points=FALSE,lines=TRUE,just=0.5)
STRIP=strip.custom(strip.names = TRUE,sep=" = ", 
                   var.name = c(expression(mu),expression(z[o])),bg=c("white"))

save_image <- function(stat,img)
{
    trellis.device(device="pdf",file=paste(stat,".pdf",sep=""),height=10,width=10)
    trellis.par.set("fontsize",list(text=18))
    print(img)
    dev.off()
}

Plot = xyplot(mLR ~ scaled_time| as.factor(mu)*as.factor(opt),#:as.factor(sweep_type),
              group=sweep_type,lwd=2,
                  type='l',data=meanLRdata,
                  par.settings=standard.theme("pdf",color=FALSE),
                  auto.key=KEY,
                  xlab="Time since optimum shift (units of N generations)",
                  ylab="Mean composite LR statistic",
                  scales=list(x=list(at=seq(-0.4,1.0,0.2),rot=45),alternating=F),
                  strip=STRIP,xlim=c(-0.5,1.0)
                  )
save_image('MeanSweeDLRLargeEffectOnly',Plot)

# nreps = rbind(hard_only,soft_only) %>% 
#     group_by(opt,mu,sweep_type) %>% 
#     summarise(n=n_distinct(repid))
# print(nreps)

# power = rbind(hard_only,soft_only) %>%
#     left_join(nreps,by=c('opt','mu','sweep_type')) %>%
#     na.omit() %>%
#     group_by(generation,opt,mu,sweep_type) %>%
#     summarise(p05=sum(sig05)/max(n),n=max(n)) %>%
#     mutate(scaled_time=(generation-5e4)/5e3)
# 
# print(unique(power$sweep_type))
# print(unique(power$n))
print(unique(meanLRdata$n))
Plot = xyplot(s05/n ~ scaled_time| as.factor(mu)*as.factor(opt),#:as.factor(sweep_type),
              group=sweep_type,lwd=2,
                  type='l',data=meanLRdata,
                  par.settings=standard.theme("pdf",color=FALSE),
                  auto.key=KEY,
                  xlab="Time since optimum shift (units of N generations)",
                  ylab="Mean composite LR statistic",
                  scales=list(x=list(at=seq(-0.4,1.0,0.2),rot=45),alternating=F),
                  strip=STRIP,xlim=c(-0.5,1)
                  )
save_image('SweeDPowerLargeEffectOnly',Plot)
