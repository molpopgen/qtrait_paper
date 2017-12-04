library(readr)
library(dplyr)
library(lattice)

sweed = read_delim("raw_sweed.txt.gz",delim=" ")
fixations = read_delim("raw_fixations.txt.gz",delim=" ")

fixations = fixations %>%
    mutate(recent = as.integer(ifelse(ftime >= 50000 & ftime < 51000,1,0))) %>%
    select(locus,repid,opt,mu,recent)

sweed_with_recent_fixations = left_join(sweed,fixations,by=c("opt","repid","locus","mu"))

sweed_with_recent_fixations = subset(sweed_with_recent_fixations,recent == 1)

x = sweed_with_recent_fixations %>% group_by(generation,opt,mu) %>%
    summarise(meanLR = mean(LR)) %>%
    mutate(scaled_time = (generation-50000)/50000)

STRIP=strip.custom(strip.names = TRUE, 
                   var.name = c(expression(z[o]),expression(mu)),bg=c("white"))
LRplot = xyplot(meanLR ~ scaled_time| as.factor(opt)*as.factor(mu),
                type='l',data=x,
                  xlab="Time since optimum shift (units of N generations)",
                  ylab="Mean composite likelihood ratio statistic",
                  scales=list(cex=1),
                  strip=STRIP)

trellis.device(device="pdf",file="SweeD_LR_recent_fixations.pdf",height=10,width=10, color=FALSE)
print(LRplot)
dev.off()
