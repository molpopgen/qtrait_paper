library(dplyr)
library(ggplot2)
library(data.table)

n=commandArgs(trailing=TRUE)

x=read.table(n[1],header=T)

##... make the plot right here.
sigEnames = list('0.1' = expression(paste(sigma["E"]," = 0.1")),
    '0.25' = expression(paste(sigma["E"]," = 0.25")))

my_labeller=function(x,y)
    {
        if(x == "r")
            {
                return(paste("r = ",as.character(y)))
            }
        return(sigEnames[as.character(y)])
    }

x.p=ggplot(x,aes(x=4*mu,y=VG)) +
    geom_point(aes(shape=as.factor(sigmu))) +
    labs(shape=expression(sigma[mu])) + 
    facet_grid(r ~ sige) + #,labeller=my_labeller) +
    geom_abline(intercept=0,slope=1,linetype="dashed") +
    xlab(expression(paste("Expected V(G)"))) +
    ylab(expression(paste("Mean simulated V(G)"))) +
    theme(legend.position="top")

ggsave(n[2])
