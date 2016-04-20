library(feather)
library(dplyr)
library(ggplot2)

x=read_feather('H2_0.2_OPT_0.5_mu_0.001.traj.rep1.feather')

                                        #Filter out sites w/sojourn times < G generations
G=5
minfreq=0.1
y=x%>%group_by(esize,pos,origin)%>%filter(length(freq)>=G & max(freq)>0.1)

rawData1=read_feather("H2_0.2_OPT_0.5_mu_0.001.popstats.feather")
rawData1SS=subset(rawData1,rep==0&(stat=="tbar"|stat=="VG"|stat=="varw"|stat=="wbar"))

#plot(y$generation,y$freq,xlim=c(9000,15000))
#lines(rawData1SS$generation[rawData1SS$stat=="tbar"],rawData1SS$value[rawData1SS$stat=="tbar"],col="red")
                                        #lines(rawData1SS$generation[rawData1SS$stat=="tbar"],rawData1SS$value[rawData1SS$stat=="wbar"],col="purple")
p = ggplot(y) +
    geom_line(aes(x=generation,y=freq,color=esize,size=origin))
