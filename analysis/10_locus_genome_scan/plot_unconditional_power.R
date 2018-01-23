library(readr)
library(lattice)
library(dplyr)
library(stringr)
library(viridis)

p = read_delim('power_tidied.txt.gz',delim=" ")

save_image <- function(stat,img)
{
    trellis.device(device="pdf",file=paste(stat,".pdf",sep=""),height=10,width=10)
    trellis.par.set("fontsize",list(text=18))
    print(img)
    dev.off()
}

KEY=list(space="top",columns=3,title="Distance from window with causal mutations.",
         cex.title=1,points=FALSE,lines=TRUE,just=0.5)
STRIP=strip.custom(strip.names = TRUE,sep=" = ", 
                   var.name = c(expression(mu),expression(z[o])),bg=c("white"))

XLIM=c(-0.5,4)
for(SI in unique(p$stat))
{
    s = p %>% filter(stat==SI,alpha==0.05) %>% mutate(dist=as.integer(abs(window-5))) %>%
        group_by(generation,mu,opt,stat,dist) %>%
        summarise(mval=mean(value)) %>%
        mutate(scaled_time = (generation-5e4)/5e3)

    COLORS=viridis(length(unique(as.factor(s$dist))))
    PLOT = xyplot(mval~scaled_time | as.factor(mu)*as.factor(opt),
                  group=dist,data=s,type='l',
                  xlim=XLIM,
                  par.settings=simpleTheme(col=COLORS),
                  auto.key=KEY,
                  scales=list(cex=1,alternating=F),
                  strip=STRIP,
                  xlab="Time since optimum shift (units of N generations)",
                  ylab="Power (averaged across unlinked loci)")


    STAT=str_to_upper(SI)
    FN=paste("Power",STAT,sep="")
    save_image(FN,PLOT)
}
