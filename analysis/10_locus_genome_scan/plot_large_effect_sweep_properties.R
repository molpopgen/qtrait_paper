library(readr)
library(dplyr)
library(viridis)
library(lattice)
library(latticeExtra)

fixations_raw = read_delim("raw_fixations.txt.gz",delim=" ")
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
#plot origin times of these events
fixations_unprocessed=fixations_raw%>%filter(sweep_type!="none",
           abs(s) >= 2*sqrt(2)*sqrt(mu),(g-5e4)/g < 0.01)
fixations_unprocessed = fixations_unprocessed %>%
    group_by(sweep_type) %>%
    mutate(ST = ifelse(sweep_type == 'hard',"Hard\nsweep","Soft\nsweep"))
print(range(fixations_unprocessed$s))
dp = bwplot(s~ST|as.factor(mu)*as.factor(opt),data=fixations_unprocessed,
                 auto.key=KEY,
                 strip=STRIP,
                  scales=list(cex=1,alternating=F),
                 panel = function(..., box.ratio) {
         panel.stripplot(..., jitter.data = T,pch=19,cex=0.5,alpha=0.2,col=viridis(1))
         panel.violin(..., col = "transparent", border = "grey60", varwidth = FALSE, 
             box.ratio = box.ratio)
    },
    ylab=expression(paste("Effect size (",gamma,")",sep="")))
save_image("LargeEffectSizeFixationsEsizeViolin",dp)
# p = xyplot(s~abs((g-5e4)/5e4)|as.factor(mu)*as.factor(opt),data=fixations_unprocessed,
#            groups=sweep_type,
#            auto.key=TRUE)
p = bwplot(abs((g-5e4)/5e4)~ST|as.factor(mu)*as.factor(opt),data=fixations_unprocessed,
           auto.key=KEY,
           scales=list(y=list(relation="free")),
           strip=STRIP,
                 panel = function(..., box.ratio) {
         panel.stripplot(..., jitter.data = T,pch=19,cex=0.5,alpha=0.2,col=viridis(1))
         panel.violin(..., col = "transparent", border = "grey60", varwidth = FALSE, 
             box.ratio = box.ratio)
    },
    ylab=expression("|Time| since optimum shift (units of N generations)"))
save_image("LargeEffectSizeFixationsTimeViolin",p)

