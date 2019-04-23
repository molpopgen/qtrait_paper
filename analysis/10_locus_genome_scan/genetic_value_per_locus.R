library(stringr)
library(readr)
library(tidyr)
library(dplyr)
library(lattice)
library(viridis)

process_gvalues <- function (filename)
{
    fs <- str_split(filename,regex('\\.'))[[1]]
    fs[2]<-str_replace(fs[2],'mu','')
    fs[4]<-str_replace(fs[4],'opt','')
    dbmu <- as.numeric(paste(fs[2],fs[3],sep='.'))
    dbopt <- as.numeric(paste(fs[4],fs[5],sep='.'))

    dbname = paste("../../mlocus_pickle/",filename,sep="")
    if(!file.exists(dbname))
    {
        stop("database not found")
    } 
    db = DBI::dbConnect(RSQLite::SQLite(),dbname)
    dbt = tbl(db,'data')

    g = collect(dbt) %>% 
        group_by(repid,generation) %>%
        mutate(rank = dense_rank(desc(abs(g)))) %>%
        group_by(generation,rank) %>%
        summarise(mg = mean(g),mag=mean(abs(g))) %>%
        mutate(mu=dbmu,opt=dbopt)
    DBI::dbDisconnect(db)
    g
}

files <- dir(path="../../mlocus_pickle/",pattern="*.genetic_values_per_locus.db")

d = NA
for (i in files)
{
    x = process_gvalues(i)
    if(is.na(d))
    {
        d=x
    }
    else
    {
        d=rbind(d,x)
    }
}
print(warnings())
d = d %>%
    mutate(scaled_time = (generation-5e4))

KEY=list(space="top",columns=3,title="",
         cex.title=1,points=FALSE,lines=TRUE,just=0.5)
STRIP=strip.custom(strip.names = TRUE, sep=" = ",
                   var.name = c(expression(mu),expression(z[o])),bg=c("white"))
COLORS=rev(viridis(length(unique(as.factor(d$rank)))))
p = xyplot((mg/opt)~scaled_time|as.factor(mu)*as.factor(opt),
           type='l',
           group=rank,
           strip=STRIP,
           #auto.key=KEY,
           par.settings=simpleTheme(col=COLORS,lwd=3),
           scales=list(cex=1,alternating=F),
           xlab="Generations since optimum shift",
           ylab=expression(paste("Contribution of locus to adaptation, ",over(bar(g), z[o]))),
           xlim=c(-0.1,0.25)*5e3,
           ylim=c(0.0,1.1),
           data=d)
trellis.device(device="pdf",file="MeanGeneticValuePerLocus.pdf",height=10,width=10)
trellis.par.set("fontsize",list(text=18))
print(p)
dev.off()
