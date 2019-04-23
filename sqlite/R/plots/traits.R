library(data.table)
library(lattice)
library(viridis)
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
data$scaled_time = (data$generation-50000)
XLIM=c(-0.02,0.05)*5e3
STRIP=strip.custom(strip.names=TRUE,var.name=expression(z[0]),
                   strip.levels=c(T,T),sep=" = ",bg=0,fg=0,style=1)

XLAB="Generations since optimum shift"
LWD=c(3,3,3)
LTY=c("solid","dashed","dotdash")
#Seriously, fuckR:
KEYTEXT = c(expression(paste(mu," = ",2.5 %*% 10^-4)), 
            expression(paste(mu," = ",1 %*% 10^-3)),        
            expression(paste(mu," = ",5 %*% 10^-3)))
KEY = list(
           space = "top", points =FALSE, lines = TRUE,
           text=KEYTEXT,
            columns = 2,cex=1)
COLORS=viridis(length(KEYTEXT))
# VGplot=xyplot(VG~scaled_time|opt,groups=mu,data=data,type='l',layout=c(1,3),
#        xlim=XLIM,
#        scales=list(x=list(tick.number=10)),
#        xlab=XLAB,
#        ylab="Genetic variance",
#        par.settings=simpleTheme(col=COLORS,lwd=3),
#        #par.settings=bwtheme,
#        strip=STRIP,
#        auto.key = KEY
#        )
VGplot=xyplot(VG~scaled_time|as.factor(opt),groups=mu,data=data,type='l',layout=c(1,3),
       xlim=XLIM,
       scales=list(cex=0.75,x=list(rot=45,tick.number=10),alternating=F),
       xlab=XLAB,
       ylab="Genetic variance",
       par.settings=simpleTheme(col=COLORS,lwd=3),
       #par.settings=bwtheme,
       strip=STRIP,
       auto.key = KEY,
       panel=function(x,y,group.number,...)
       {
           EVG = 4*unique(data$mu)[group.number]
           panel.xyplot(x,y,col=COLORS[group.number],lwd=3,...)
           panel.abline(h=EVG,lty="dotdash",col=COLORS[group.number],lwd=2)
       }
       )
trellis.device(device="pdf",file="slocusVG.pdf")
trellis.par.set("fontsize",list(text=18))
print(VGplot)
dev.off()
trellis.device(device="tiff",file="slocusVG.tiff")
trellis.par.set("fontsize",list(text=18))
print(VGplot)
dev.off()

zbar=xyplot(tbar~scaled_time|as.factor(opt),group=mu,data=data,type='l',layout=c(1,3),
       xlim=XLIM,
       scales=list(cex=0.75,x=list(rot=45,tick.number=10),alternating=F),
       xlab=XLAB,
       ylab=expression(bar(z)),
       par.settings=simpleTheme(col=COLORS,lwd=3),
       #par.settings=bwtheme,
       strip=STRIP,
       auto.key = KEY
       )
trellis.device(device="pdf",file="slocusZbar.pdf")

trellis.par.set("fontsize",list(text=18))
print(zbar)
dev.off()
