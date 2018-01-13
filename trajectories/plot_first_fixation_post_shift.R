library(readr)
library(lattice)

data = read_delim("first_fixation_properties.txt",delim=" ")

data$scaled_time = (data$generation - 5e4)/5e3
data$median_first_ftime = (data$median_first_ftime-5e4)/5e3
data$median50= (data$mean50-5e4)/5e3
ds = subset(data,opt==0.5)
XLAB="Time since optimum shift (units of N generations)"
STRIP=strip.custom(strip.names=TRUE,var.name=c(expression(mu),expression(z[0])),
                   strip.levels=c(T,T),sep=" = ",bg=0,fg=0,style=1)
p = xyplot(mz~scaled_time|as.factor(mu)*as.factor(opt),data=ds,
           type='l',scales =
       list(x =
            list(relation = "free")),
           xlab=XLAB,
           ylab="Mean trait value",
           strip=STRIP,
           panel=function(x, y,...){
            panel.xyplot(x,y, type='l')
               panel_mu = unique(ds$mu)[panel.number()]
               median_first_ftime = ds$median_first_ftime[ds$mu == panel_mu][1]
               m10 = ds$median10 = ds$median10[ds$mu==panel_mu][1]
               m50 = ds$median50 = ds$median50[ds$mu==panel_mu][5]
            panel.abline(v=median_first_ftime)
            panel.abline(v=m50,lty="dotted")
})


