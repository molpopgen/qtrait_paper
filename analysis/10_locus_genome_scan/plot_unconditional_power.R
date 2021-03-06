# This script is replaced with the _V2 version,
# which is mo' correct.
library(readr)
library(lattice)
library(dplyr)
library(stringr)
library(viridis)

p = read_delim('nsig_per_window_tidied.txt.gz',delim=" ")

save_image <- function(stat,img)
{
    trellis.device(device="pdf",file=paste(stat,".pdf",sep=""),height=10,width=10)
    trellis.par.set("fontsize",list(text=18))
    print(img)
    dev.off()
}

# KEY=list(space="top",columns=3,title="Distance from window with causal mutations.",
#          cex.title=1,points=FALSE,lines=TRUE,just=0.5)
STRIP=strip.custom(strip.names = TRUE,sep=" = ", 
                   var.name = c(expression(mu),expression(z[o])),bg=c("white"))

XLIM=c(-0.5,4)
for(SI in unique(p$stat))
{
    s = p %>% filter(stat==SI,alpha==0.05) %>% mutate(dist=as.integer(abs(window-5))) %>%
        group_by(generation,mu,opt,stat,dist) %>%
        summarise(total=sum(value)) %>%
        mutate(scaled_time = (generation-5e4)/5e3,
               # The adjustment is due to all windows
               # other than the selected one have
               # stats that are pooled via the dist
               # factor, meaning 2x the number
               # of observations.
               adjustment = ifelse(dist == 0,1,2))

    ICOLORS=viridis(length(unique(as.factor(s$dist))))
    COLORS=array()
    COLORS[1] = rgb(as.integer(col2rgb(ICOLORS[1])[1,])/255,
                as.integer(col2rgb(ICOLORS[1])[2,])/255,
                as.integer(col2rgb(ICOLORS[1])[3,])/255,
                alpha=1/length(ICOLORS))
    for(i in 2:length(ICOLORS))
    {
        ncolor = rgb(as.integer(col2rgb(ICOLORS[i])[1,])/255,
                    as.integer(col2rgb(ICOLORS[i])[2,])/255,
                    as.integer(col2rgb(ICOLORS[i])[3,])/255,
                    alpha=i/length(ICOLORS))
        COLORS[i]=ncolor
    }
    nlines = length(unique(s$dist))
    KEY=list(space="top",columns=3,
             title="Distance from window with causal mutations.",
             cex.title=1,#points=FALSE,
             lines=list(lwd=rep(3,length(COLORS)),col=rev(COLORS)),
             just=0.5,
             text=list(as.character(sort(unique(s$dist)))))

    PLOT = xyplot(total/(adjustment*256*11)~scaled_time | as.factor(mu)*as.factor(opt),
                  lwd=3,data=s,type='l',
                  group=factor(dist,levels=rev(sort(unique(s$dist)))),
                  #group=dist,data=s,type='l',
                  xlim=XLIM,
                  par.settings=simpleTheme(col=COLORS),
                  key=KEY,
                  scales=list(cex=1,alternating=F),
                  strip=STRIP,
                  xlab="Time since optimum shift (units of N generations)",
                  ylab=expression(paste("Power at ",alpha," = 0.05")))


    STAT=str_to_upper(SI)
    FN=paste("Power",STAT,sep="")
    save_image(FN,PLOT)
}

