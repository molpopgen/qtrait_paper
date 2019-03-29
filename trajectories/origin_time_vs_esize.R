library(stringr)
library(dplyr)
library(lattice)

process <- function(f)
{
    s=str_split(f,'_')
    opt=as.numeric(s[[1]][4])
    mu=as.numeric(s[[1]][6])

    print(f)
    db = DBI::dbConnect(RSQLite::SQLite(),f)
    dbt = tbl(db,'freqs')

    q = dbt %>% filter(origin < 60000)

    fixations = collect(dbt) %>% 
        group_by(repid,origin,pos,esize) %>%
        summarise(n=n()) %>% # this is the sojourn time
        # group_by(origin) %>%
        # summarise(mesize=mean(esize),mean_ftime = mean(n-origin)) %>%
        mutate(opt=opt,mu=mu)


    print(fixations)

    DBI::dbDisconnect(db)

    fixations
}


files = dir('.',"_traj_merged_recent_fixations.db$")

data = NA
I=0
for(i in files)
{
    x = process(i)
    if(I==0)
    {
        data=x
    }
    else
    {
        data = rbind(data,x)
    }
    I=I+1
}
KEY=list(space="top",columns=3,title="",
         cex.title=1,points=FALSE,lines=TRUE,just=0.5)
STRIP=strip.custom(strip.names = TRUE,sep=" = ",
                   var.name = c(expression(mu),expression(z[o])),bg=c("white"))
#COLORS=c("black")

data$scaled_origin = (data$origin - 5e4)# /5e3

p = xyplot(esize~scaled_origin|as.factor(mu)*as.factor(opt),data=data,
            mu = data$mu,sojourn=data$n,optima=data$opt,xlim=c(-1250,1250),
            par.settings=simpleTheme(cex=1,alpha.points=0.5),
            auto.key=KEY, xlab="Origin time of fixation (generations since optimum shift)",
            ylab=expression(paste("Effect size (",gamma,")")),
            scales=list(cex=1,alternating=F,x=list(rot=45)),
            strip=STRIP,
            panel=function(x, y,mu,sojourn,optima,...,subscripts){
            panel_sojourn = sojourn[subscripts]
            gamma_overshoots = unique(optima[subscripts])/2
            panel.xyplot(x[x<0 & x+panel_sojourn<0 ],y[x<0 & x+panel_sojourn<0],col="black",pch=22)
            panel.xyplot(x[x<0 & x+panel_sojourn>0 ],y[x<0 & x+panel_sojourn>0],col="red",pch=16)
            panel.xyplot(x[x>0],y[x>0],pch=15)#,col="blue")
            panel_mu = unique(mu[subscripts])
            # ghat = 2*sqrt(2)*sqrt(panel_mu)
            # for (gcrit in c(10, 100, 1000)){
            #     big = sqrt(gcrit/(5000))
            #     panel.abline(h=big,lwd=2,lty='dotdash')
            #     panel.abline(h=-big,lwd=2,lty='dotdash')
            # }
            panel.abline(v=0.0,lty="dotted",lwd=2)
})



trellis.device(device="pdf",file="OriginTimeVsEsizeFixations.pdf",height=10,width=10)
trellis.par.set("fontsize",list(text=18))
print(p)
dev.off()
