library(dplyr)
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
        dbname = paste(outer[i],"_",inner[j],"_qtrait.sqlite3",sep="")
        db = dbConnect(RSQLite::SQLite(), dbname)
        t = tbl(db, "data")

        df = collect(t) %>% group_by(generation) %>%
            summarise(mz = mean(zbar), vz = var(zbar)) %>%
            mutate(mu=mu[j], rho=rho[i])

        data = bind_rows(data, df)

        dbDisconnect(db)
    }
}

data = data %>%
    mutate(tsince = generation - 5e4)

COLORS=viridis(3)
KEY=list(space="top",columns=3,
         title=expression(paste("Scaled recombination rate, ",rho)),
         cex.title=0.75,#points=FALSE,
         lines=list(lwd=rep(3,length(COLORS)),col=COLORS),
         just=0.5,
         text=list(as.character(sort(rho))))
STRIP=strip.custom(strip.names = TRUE,sep=" = ", 
                   var.name = expression(mu), bg=c("white"))
plt = xyplot(mz~tsince|as.factor(mu), data=data, groups=as.factor(rho), type='l',
             par.settings=simpleTheme(col=COLORS),
             lwd=3, xlim=c(-10,110),
             scales=list(cex=0.75,alternating=F,x=list(rot=45,at=seq(0,100,20))),
             layout=c(3,1),
             key=KEY, strip=STRIP,
             xlab="Generations since optimum shift",
             ylab="Mean trait value")
trellis.device(device="pdf",file="MeanGeneticValueVaryRho.pdf",height=10,width=10)
trellis.par.set("fontsize",list(text=18))
print(plt)
dev.off()
plt = xyplot(vz~tsince|as.factor(mu), data=data, groups=as.factor(rho), type='l',
             par.settings=simpleTheme(col=COLORS),
             lwd=3, xlim=c(-10,110),
             scales=list(cex=0.75,alternating=F,x=list(rot=45,at=seq(0,100,20))),
             layout=c(3,1),
             key=KEY, strip=STRIP,
             xlab="Generations since optimum shift",
             ylab="Mean genetic variance")
trellis.device(device="pdf",file="MeanVGVaryRho.pdf",height=10,width=10)
trellis.par.set("fontsize",list(text=18))
print(plt)
dev.off()
