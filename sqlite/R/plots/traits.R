library(data.table)
library(lattice)
library(ggplot2)

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
bwtheme <- standard.theme("pdf", color=FALSE)

pdf("slocusVG.pdf")
print(paste("z","=",as.numeric(optvals)))
data$scaled_time = (data$generation-50000)/5000.0
xyplot(VG~scaled_time|as.factor(opt),group=mu,data=data,type='l',layout=c(1,3),
       xlim=c(-0.1,0.25),
       scales=list(x=list(tick.number=10)),
       xlab="Time since optimum shift (units of N generation)",
       ylab=expression(V[G]),
       par.settings=bwtheme,
       strip=strip.custom(var.name=expression(z[0]),strip.levels=c(T,T),sep=" = ",bg=0),
       lwd=c(3,3,3)
       )
dev.off()
