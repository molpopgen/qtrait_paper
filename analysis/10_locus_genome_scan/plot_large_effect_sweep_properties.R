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
#plot origin times of these events
fixations_unprocessed=fixations_raw%>%filter(sweep_type!="none",
           abs(s) >= 2*sqrt(2)*sqrt(mu),(g-5e4)/g < 0.01)
dp = densityplot(~s|as.factor(mu)*as.factor(opt),data=fixations_unprocessed,groups=sweep_type,
                 auto.key=TRUE)
save_image("esizes_large_effect_sweeps",dp)
p = xyplot(s~abs((g-5e4)/5e4)|as.factor(mu)*as.factor(opt),data=fixations_unprocessed,
           groups=sweep_type,
           auto.key=TRUE)
save_image("esizes_vs_distance_from_shift_large_effect_sweeps",p)

