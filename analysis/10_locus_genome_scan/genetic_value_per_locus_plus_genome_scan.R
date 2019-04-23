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
        mutate(locus=locus,rank = dense_rank(desc(abs(g)))) %>%
        # group_by(generation,rank,locus) %>%
        # summarise(mg = mean(g),mag=mean(abs(g))) %>%
        mutate(mu=dbmu,opt=dbopt)
    DBI::dbDisconnect(db)
    g
}

proess_gscan <- function(filename)
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
        filter(window == 5) %>%  #window w/selected sites
        # group_by(locus,generation) %>%
        # summarise(meanD=mean(tajd),meanH=mean(hprime)) %>%
        select(repid,locus,generation,tajd,hprime) %>%
        mutate(opt=dbopt,mu=dbmu)
    DBI::dbDisconnect(db)
    g
}

files <- dir(path="../../mlocus_pickle/",pattern="*.genetic_values_per_locus.db")

d = NA
for (i in files) #[1])
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
gscan_files <- dir(path="../../mlocus_pickle/",pattern="*.genome_scan.db")
g = NA
for (i in gscan_files) #[1])
{
    x = proess_gscan(i)
    if(is.na(g))
    {
        g=x
    }
    else
    {
        g=rbind(g,x)
    }
}

print(head(d))
print(head(g))

d = d %>% full_join(g,by=c('repid','locus','generation','opt','mu')) %>%
    group_by(generation,rank,opt,mu) %>%
    summarise(meanD=mean(tajd),meanH=mean(hprime)) %>%
    mutate(scaled_time = (generation-5e4))

# print(head(d))

#q("no")
#print(warnings())
# d = d %>%

KEY=list(space="top",columns=3,title="",
         cex.title=1,points=FALSE,lines=TRUE,just=0.5)
STRIP=strip.custom(strip.names = TRUE, sep=" = ",
                   var.name = c(expression(mu),expression(z[o])),bg=c("white"))
COLORS=rev(viridis(length(unique(as.factor(d$rank)))))
p = xyplot(meanD~scaled_time|as.factor(mu)*as.factor(opt),
           type='l',
           group=rank,
           strip=STRIP,
           #auto.key=KEY,
           par.settings=simpleTheme(col=COLORS,lwd=3),
           scales=list(cex=1,alternating=F),
           xlab="Time since optimum shift (units of N generations)",
           ylab="Mean Tajima's D",
           xlim=c(-0.1,3)*5e3,
           data=d)
trellis.device(device="pdf",file="MeanTajDPerRankedLocus.pdf",height=10,width=10)
trellis.par.set("fontsize",list(text=18))
print(p)
dev.off()

p = xyplot(meanH~scaled_time|as.factor(mu)*as.factor(opt),
           type='l',
           group=rank,
           strip=STRIP,
           #auto.key=KEY,
           par.settings=simpleTheme(col=COLORS,lwd=3),
           scales=list(cex=1,alternating=F),
           xlab="Time since optimum shift (units of N generations)",
           ylab="Mean H'",
           xlim=c(-0.1,3)*5e3,
           data=d)
trellis.device(device="pdf",file="MeanHprimePerRankedLocus.pdf",height=10,width=10)
trellis.par.set("fontsize",list(text=18))
print(p)
dev.off()

