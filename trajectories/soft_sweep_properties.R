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

    softq = dbt %>% filter(origin < 50000, generation >= 50000) 

    soft = collect(softq) %>% 
        group_by(repid,origin,pos,esize) %>%
        summarise(n=n(),freq_at_shift = unique(freq[generation == 50000])) %>% # this is the sojourn time
        mutate(opt=opt,mu=mu)

    DBI::dbDisconnect(db)

    soft
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
STRIP=strip.custom(strip.names = TRUE, 
                   var.name = c(expression(z[o]),expression(mu)),bg=c("white"))
#COLORS=c("black")

p = xyplot(esize~freq_at_shift|as.factor(mu)*as.factor(opt),data=data,
           mu = data$mu,
          par.settings=simpleTheme(pch=16,cex=0.5,alpha.points=0.5),
          auto.key=KEY, xlab="Mutation frequency at optimum shift",
          ylab=expression(paste("Effect size (",gamma,")")),
          scales=list(cex=1,alternating=F,x=list(rot=45)),
          strip=STRIP,
           panel=function(x, y,mu,...,subscripts){
            panel.xyplot(x,y)
            panel_mu = unique(mu[subscripts])
           ghat = 2*sqrt(2)*sqrt(panel_mu)
            panel.abline(h=ghat)
    })



trellis.device(device="pdf",file="SoftSweepFreqEsize.pdf",height=10,width=10)
print(p)
dev.off()
