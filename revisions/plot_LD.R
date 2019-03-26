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
        dbname = paste(outer[i],"_",inner[j],"_ld.sqlite3",sep="")
        db = dbConnect(RSQLite::SQLite(), dbname)
        t = tbl(db, "data")

        df = collect(t) %>% group_by(generation) %>%
            summarise(mintraD=mean(intralocus_D),
                      minterD=mean(interlocus_D)) %>%
            mutate(mu=mu[j], rho=rho[i])

        data = bind_rows(data, df)

        dbDisconnect(db)
    }
}

COLORS=viridis(3)
data = data %>% filter(mu > 2.5e-4) %>%
    mutate(tsince = generation - 5e4) %>%
    mutate(color = NA)

data$color[data$rho == 1e2] <- COLORS[1]
data$color[data$rho == 1e3] <- COLORS[2]
data$color[data$rho == 1e4] <- COLORS[3]

print(COLORS)
KEY=list(space="top",columns=3,
         title=expression(paste("Scaled recombination rate, ",rho)),
         cex.title=0.75,#points=FALSE,
         lines=list(lwd=rep(3,length(COLORS)),col=COLORS),
         just=0.5,
         text=list(as.character(sort(rho))))
STRIP=strip.custom(strip.names = TRUE,sep=" = ", 
                   var.name = expression(mu), bg=c("white"))
plt = xyplot(mintraD ~ tsince|as.factor(mu), data=data, 
             minterD=data$minterD,
             color=data$color,
             rho=data$rho,
            group=as.factor(rho),
            key=KEY, strip=STRIP,
            xlab="Generations since optimum shift",
            ylab="Mean disequilibrium statistic, D",
             scales=list(cex=0.75,alternating=F,x=list(rot=45,at=seq(0,200,50))),
            panel=function(x,y,minterD,color,rho,...,subscripts)
            {
                rs = rho[subscripts]
                cs = color[subscripts]
                mDs = minterD[subscripts]
                for (ur in unique(rs))
                {
                    theseindexes = which(rs == ur)
                    thiscolor = unique(cs[theseindexes])
                    panel.xyplot(x[theseindexes],y[theseindexes],type='l',lwd=3,
                                 col.line=thiscolor)
                    panel.xyplot(x[theseindexes],minterD[theseindexes],type='l',lwd=3,
                                 lty='dotdash',col.line=thiscolor)
                }
            })


trellis.device(device="pdf",file="MeanLDVaryRho.pdf",height=10,width=10)
trellis.par.set("fontsize",list(text=18))
print(plt)
dev.off()
