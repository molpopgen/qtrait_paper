library(dplyr)
library(stringr)
library(lattice)
library(readr)

x = read_delim("raw_sweed.txt.gz",delim=" ")

x = x %>% group_by(generation,mu,opt) %>%
    summarise(meanLR = mean(LR), mean_alpha = mean(alpha)) %>%
    mutate(scaled_time = (generation - 50000))

STRIP=strip.custom(strip.names = TRUE,sep=" = ", 
                   var.name = c(expression(mu),expression(z[o])),bg=c("white"))
LRplot = xyplot(meanLR ~ scaled_time| as.factor(mu)*as.factor(opt),
                type='l',data=x,
                  par.settings=simpleTheme(),
                  xlab="Generations since optimum shift",
                  ylab="Mean composite likelihood ratio statistic",
                  xlim=c(-0.5,1.0)*5e3,
                  scales=list(x=list(at=seq(-0.4,0.8,0.2)*5e3,rot=45),alternating=F),
                  strip=STRIP)

trellis.device(device="pdf",file="SweeDLRRaw.pdf",height=10,width=10, color=FALSE)
trellis.par.set("fontsize",list(text=18))
print(LRplot)
dev.off()

alpha_plot = xyplot(mean_alpha ~ scaled_time| as.factor(mu)*as.factor(opt),
                type='l',data=x,
                  xlab="Generations since optimum shift",
                  ylab="Mean estimate of 2Ns",
                  scales=list(cex=1),
                  strip=STRIP)

trellis.device(device="pdf",file="SweedAlphaRaw.pdf",height=10,width=10, color=FALSE)
trellis.par.set("fontsize",list(text=18))
print(alpha_plot)
dev.off()
