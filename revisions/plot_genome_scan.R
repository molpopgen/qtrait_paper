library(dplyr)
library(tidyr)
library(DBI)
library(lattice)
library(viridis)

outer = c("main_results", "lowrec", "hirec")
rho = c(1e3,1e2,1e4)
inner = c("lomu", "midmu", "himu")
mu = c(2.5e-4, 1e-3, 5e-3)

data = tibble()
for (i in 1:length(outer))
{
    for(j in 1:length(inner))
    {
        dbname = paste(outer[i],"_",inner[j],"_gscan.sqlite3",sep="")
        db = dbConnect(RSQLite::SQLite(), dbname)
        t = tbl(db, "data")

        df = collect(t) %>% 
            mutate(dist = abs(window-5)) %>%
            group_by(generation, dist) %>%
            summarise(mpi=mean(pi), mD=mean(tajd),mHp=mean(hprime)) %>%
            mutate(mu=mu[j], rho=rho[i])

        data = bind_rows(data, df)

        dbDisconnect(db)
    }
}
data = data %>% mutate(tsince = generation-5e4)


ICOLORS=viridis(length(unique(as.factor(data$dist))))
COLORS=array()
COLORS[1] = rgb(as.integer(col2rgb(ICOLORS[1])[1,])/255,
            as.integer(col2rgb(ICOLORS[1])[2,])/255,
            as.integer(col2rgb(ICOLORS[1])[3,])/255,
            alpha=1/length(ICOLORS))
for(i in 2:length(ICOLORS))
{
    ncolor = rgb(as.integer(col2rgb(ICOLORS[i])[1,])/255,
                as.integer(col2rgb(ICOLORS[i])[2,])/255,
                as.integer(col2rgb(ICOLORS[i])[3,])/255,
                alpha=i/length(ICOLORS))
    COLORS[i]=ncolor
}

KEY=list(space="top",columns=3,
         title="Distance from window with causal mutations.",
         cex.title=1,#points=FALSE,
         lines=list(lwd=rep(3,length(COLORS)),col=rev(COLORS)),
         just=0.5,
         text=list(as.character(sort(unique(data$dist)))))
STRIP=strip.custom(strip.names = TRUE,sep=" = ", 
                   var.name = c(expression(mu),expression(rho)), bg=c("white"))
plt <- xyplot(mD ~ tsince|as.factor(mu)*as.factor(rho),data=data,
              group=factor(dist,levels=rev(sort(unique(data$dist)))),
              par.settings=simpleTheme(col=COLORS),
              lwd=3,
              type='l',
              key=KEY,
              scales=list(cex=1,alternating=F,x=list(rot=45,at=seq(0,200,50))),
              strip=STRIP,
              xlab="Generations since optimum shift",
              ylab="Mean Tajima's D")

trellis.device(device='pdf',file='MeanTajimasDVaryRho.pdf',height=10,width=10)
trellis.par.set("fontsize",list(text=18))
print(plt)
dev.off()

plt <- xyplot(mHp ~ tsince|as.factor(mu)*as.factor(rho),data=data,
              group=factor(dist,levels=rev(sort(unique(data$dist)))),
              par.settings=simpleTheme(col=COLORS),
              lwd=3,
              type='l',
              key=KEY,
              scales=list(cex=1,alternating=F,x=list(rot=45,at=seq(0,200,50))),
              strip=STRIP,
              xlab="Generations since optimum shift",
              ylab="Mean H'")

trellis.device(device='pdf',file='MeanHprimeVaryRho.pdf',height=10,width=10)
trellis.par.set("fontsize",list(text=18))
print(plt)
dev.off()
