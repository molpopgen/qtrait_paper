library(data.table)
library(lattice)
library(ggplot2)
source("latticebw.R")

muvals = c("0.00025", "0.001", "0.005")
optvals = c("0.1","0.5","1")

data = NA
for (m in muvals)
{
    for(o in optvals)
    {
        fn = paste("../H2_1.0_OPT_",o,"_mu_",m,"_sigmu0.25_stats.gz",sep="")
        x = fread(paste("gunzip -c ",fn,sep=""))
        if( is.na(data) )
        {
            data = x
        }
        else
        {
            data = rbind(data,x)
        }
    }
}

#Lattice graphics, FTW
data$scaled_time = (data$generation-50000)/5000.0
XLIM=c(-0.02,0.1)
STRIP=strip.custom(var.name=expression(z[0]),strip.levels=c(T,T),sep=" = ",bg=0,fg=0,style=1)
XLAB="Time since optimum shift (units of N generation)"
LWD=c(3,3,3)
LTY=c("solid","dashed","dotdash")
KEYTEXT = c(expression(paste(mu," = ",2.5e-4)),
            expression(paste(mu," = ",1e-3)),        
            expression(paste(mu," = ",5e-3)))
KEY = list(
           space = "top", points = FALSE, lines = TRUE,
           text=KEYTEXT,
            columns = 3)
VGplot=xyplot(VG~scaled_time|opt,group=mu,data=data,type='l',layout=c(1,3),
       xlim=XLIM,
       scales=list(x=list(tick.number=10)),
       xlab=XLAB,
       ylab="Genetic variance",
       par.settings=bwtheme,
       strip=STRIP,
       auto.key = KEY
       )
pdf("slocusVG.pdf")
VGplot
dev.off()
tiff("slocusVG.tiff")
VGplot
dev.off()

pdf("slocusZbar.pdf")
xyplot(tbar~scaled_time|opt,group=mu,data=data,type='l',layout=c(1,3),
       xlim=XLIM,
       scales=list(x=list(tick.number=10)),
       xlab=XLAB,
       ylab=expression(bar(z)),
       par.settings=bwtheme,
       strip=STRIP,
       auto.key = KEY,
       )
dev.off()
