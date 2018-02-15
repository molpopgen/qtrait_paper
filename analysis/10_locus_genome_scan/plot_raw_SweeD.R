library(dplyr)
library(stringr)
library(lattice)
library(readr)

x = read_delim("raw_sweed.txt.gz",delim=" ")

x = x %>% group_by(generation,mu,opt) %>%
    summarise(meanLR = mean(LR), mean_alpha = mean(alpha)) %>%
    mutate(scaled_time = (generation - 50000)/50000)

STRIP=strip.custom(strip.names = TRUE, 
                   var.name = c(expression(mu),expression(opt)),bg=c("white"))
LRplot = xyplot(meanLR ~ scaled_time| as.factor(mu)*as.factor(opt),
                type='l',data=x,
                  xlab="Time since optimum shift (units of N generations)",
                  ylab="Mean composite likelihood ratio statistic",
                  scales=list(cex=1),
                  strip=STRIP)

trellis.device(device="pdf",file="SweeDLRRaw.pdf",height=10,width=10, color=FALSE)
print(LRplot)
dev.off()

alpha_plot = xyplot(mean_alpha ~ scaled_time| as.factor(mu)*as.factor(opt),
                type='l',data=x,
                  xlab="Time since optimum shift (units of N generations)",
                  ylab="Mean estimate of 2Ns",
                  scales=list(cex=1),
                  strip=STRIP)

trellis.device(device="pdf",file="SweedAlphaRaw.pdf",height=10,width=10, color=FALSE)
print(alpha_plot)
dev.off()
