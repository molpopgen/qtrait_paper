#Plots for bio of genomes poster
library(feather)
library(ggplot2)

x=read_feather('popstats.feather')

print(names(x))
print(head(x))

x1=subset(x,stat=="VG"|stat=="max_expl")#|stat=="leading_e"|stat=="leading_q")

y=subset(x1,H2=="0.2"&mu==0.001)

##Means, etc., of "quant-gen" params
posterFigure1 = ggplot(y) + geom_line(aes(x=(generation-10000)/1e3,y=value,color=stat),size=1 )+
    facet_wrap(~opt,nrow=2,scales="free") +
        xlab("Time since optimum shift (units of N generations)") +
            ylab("Mean value of statistic") +
                scale_color_discrete(labels=c(expression(paste("max[2p(1-p)",e^2,"]")),expression(VG))) +
                theme_bw(base_size=14) +
                    theme(legend.position="top") +
                        coord_cartesian(xlim=c(-1,1.5))
ggsave('posterFigure1.pdf',posterFigure1)
                                        #x1=subset(x,stat=="VG"|stat=="max_expl"|stat=="tbar"|stat=="leading_e"|stat=="leading_q"|stat=="wbar")
x1=subset(x,stat=="VG"|stat=="varw"|stat=="tbar"|stat=="wbar")
x1$value[x1$stat=="tbar"]=x1$value[x1$stat=="tbar"]*2
x1$value[x1$stat=="VG"]=x1$value[x1$stat=="VG"]*10
x1$value[x1$stat=="varw"]=x1$value[x1$stat=="varw"]*200
y=subset(x1,opt==0.5&mu==0.001)

##Pop'n reaches optimum, but there is still variation in fitness, etc.
posterFigure2 = ggplot(y) +
    geom_line(aes(x=(generation-10000)/1e3,y=value,color=stat),size=2) +
        #facet_wrap(~opt) +
            coord_cartesian(xlim=c(-0.5,1.5)) +
                xlab("Time since optimum shift (units of N generations)") +
                    ylab("Mean value of statistic") +
                        geom_vline(xintercept=0.25) +
                            geom_segment(aes(xend=0.26,x=0.5,y=0.75,yend=0.75),arrow=arrow(length=unit(0.5,"cm"))) +
                                annotate("text",x=1.0,y=0.75,label="Mean time to reach\noptimum trait value") +
                                    theme_bw(base_size=14) +
                                        scale_color_discrete(labels=c("2*[Mean trait value]",
                                                                 "200*Var(fitness)",
                                                                 "10*VG",
                                                                 expression(bar(w))))
ggsave('posterFigure2.pdf',posterFigure2)


#This is a plot of the "quant-gen" parameters for replicate 0
rawData1=read_feather("H2_0.2_OPT_0.5_mu_0.001.popstats.feather")
rawData1SS=subset(rawData1,rep==0&(stat=="tbar"|stat=="VG"|stat=="varw"|stat=="wbar"))
rawData1SS$value[rawData1SS$stat=="tbar"]=2*rawData1SS$value[rawData1SS$stat=="tbar"]
rawData1SS$value[rawData1SS$stat=="VG"]=10*rawData1SS$value[rawData1SS$stat=="VG"]
rawData1SS$value[rawData1SS$stat=="varw"]=200*rawData1SS$value[rawData1SS$stat=="varw"]
rawData1SSplot = ggplot(rawData1SS) +
    geom_line(aes(x=(generation-10000)/1e3,y=value,color=stat),size=1.25) +
        coord_cartesian(xlim=c(-0.5,1.5)) +
            xlab("Time since optimum shift (units of N generations)") +
                ylab("Mean value of statistic") +
                    geom_vline(xintercept=0.25) +
                        theme_bw(base_size=14) +
                            scale_color_discrete(labels=c("2*[Mean trait value]",
                                                     "200*Var(fitness)",
                                                     "10*VG",
                                                     expression(bar(w))))
ggsave("posterFigure3.pdf",rawData1SSplot)
