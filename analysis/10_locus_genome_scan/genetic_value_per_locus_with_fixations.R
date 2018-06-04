library(stringr)
library(readr)
library(tidyr)
library(dplyr)
library(lattice)
library(viridis)

process_gvalues <- function (filename, fixations)
{
    fs <- str_split(filename,regex('\\.'))[[1]]
    fs[2]<-str_replace(fs[2],'mu','')
    fs[4]<-str_replace(fs[4],'opt','')
    dbmu <- as.numeric(paste(fs[2],fs[3],sep='.'))
    dbopt <- as.numeric(paste(fs[4],fs[5],sep='.'))

    f = fixations %>% filter(opt==dbopt,mu==dbmu)
    dbname = paste("../../mlocus_pickle/",filename,sep="")
    if(!file.exists(dbname))
    {
        stop("database not found")
    } 
    db = DBI::dbConnect(RSQLite::SQLite(),dbname)
    dbt = tbl(db,'data')

    g = collect(dbt) %>% mutate(opt=dbopt,mu=dbmu) %>%
        group_by(repid,generation) %>% mutate(rank = dense_rank(desc(abs(g)))) %>%
        inner_join(f,by=c('generation','repid','locus','mu','opt')) %>%
        group_by(generation,rank,sweep_type) %>%
        summarise(meanD = mean(tajd),meanH=mean(hprime),mg=mean(g)) %>%
        mutate(opt=dbopt,mu=dbmu) %>%
        filter(sweep_type != 'none')
    # print(paste(unique(g$opt),unique(g$mu)))
    # print(g)
    # print(names(g))
    # print("leaving")
    # q("no")
    # g = collect(dbt) %>%
    #     group_by(repid,generation) %>%
    #     mutate(rank = dense_rank(desc(abs(g)))) %>%
    #     mutate(opt=dbopt,mu=dbmu) %>%
    #     full_join(f,by=c('generation','opt','mu','repid','locus')) %>%
    #     drop_na() %>%
    #     group_by(generation,rank,sweep_type) %>%
    #     summarise(meanD = mean(tajd),meanH=mean(hprime),mg=mean(g)) %>%
    #     mutate(opt=dbopt,mu=dbmu)
        
    DBI::dbDisconnect(db)
    g
}

fixations_db=DBI::dbConnect(RSQLite::SQLite(),"genome_scan_with_large_effect_sweep_counts.sqlite3")
fixations_t = tbl(fixations_db,'data')
fixations = collect(fixations_t)

fixations = fixations %>%
    mutate(sweep_type = ifelse(nhard > 0 & nsoft == 0,
                               'New mutation',
                               ifelse(nsoft>0&nhard==0,
                                      'Standing var.','none'))) %>%
    # Retain genome scan info 
    # only for window w/selected
    # fixations
    filter(window==5)

DBI::dbDisconnect(fixations_db)

files <- dir(path="../../mlocus_pickle/",pattern="*.genetic_values_per_locus.db")

d = NA
for (i in files) #[1:1])
{
    print(paste("processing",i))
    x = process_gvalues(i, fixations)
    if(is.na(d))
    {
        d=x
    }
    else
    {
        d=rbind(d,x)
    }
}
d = d %>%
    mutate(scaled_time = (generation-5e4)/5e3)

KEY=list(space="top",columns=3,title="",
         cex.title=1,points=FALSE,lines=TRUE,just=0.5)
STRIP=strip.custom(strip.names = TRUE, sep=" = ",
                   var.name = c(expression(mu),expression(z[o])),bg=c("white"))
COLORS=rev(viridis(length(unique(as.factor(d$rank)))))
p = xyplot((mg/opt)~scaled_time|as.factor(mu)*as.factor(opt):as.factor(sweep_type),
           type='l',
           group=rank,
           strip=STRIP,
           #auto.key=KEY,
           par.settings=simpleTheme(col=COLORS,lwd=3),
           scales=list(cex=0.75,alternating=F),
           xlab="Time since optimum shift (units of N generations)",
           ylab=expression(paste("Contribution of locus to adaptation, ",over(bar(g), z[o]))),
           xlim=c(-0.1,0.25),
           ylim=c(0.0,1.1),
           data=subset(d,opt > 0.1)
           )
trellis.device(device="pdf",file="MeanGeneticValuePerLocusLargeEffect.pdf",height=10,width=10)
trellis.par.set("fontsize",list(text=18))
print(p)
dev.off()
d2 = d %>% filter(opt>0.1,rank<3) %>% mutate(G = ifelse(sweep_type == 'New mutation',ifelse(rank==1,"N1","N2"),ifelse(rank==1,"S1","S2")))
p = xyplot(meanD~scaled_time|as.factor(mu)*as.factor(opt),#:as.factor(sweep_type),
           type='l',
           group=G,
           strip=STRIP,
           par.settings=simpleTheme(col=c(COLORS[1],COLORS[2],COLORS[1],COLORS[2]),
                                    lty=c("solid","solid","dashed","dashed"),lwd=3),
           scales=list(cex=0.75,alternating=F),
           xlab="Time since optimum shift (units of N generations)",
           ylab="Mean Tajima's D",
           xlim=c(-0.1,3),
           #data=subset(d,opt > 0.1&rank<3)
           data=d2
           )
trellis.device(device="pdf",file="MeanTajDPerLocusLargeEffect.pdf",height=10,width=10)
trellis.par.set("fontsize",list(text=18))
print(p)
dev.off()
p = xyplot(meanH~scaled_time|as.factor(mu)*as.factor(opt),#:as.factor(sweep_type),
           type='l',
           group=G,
           strip=STRIP,
           #auto.key=KEY,
           par.settings=simpleTheme(col=c(COLORS[1],COLORS[2],COLORS[1],COLORS[2]),
                                    lty=c("solid","solid","dashed","dashed"),lwd=3),
           scales=list(cex=0.75,alternating=F),
           xlab="Time since optimum shift (units of N generations)",
           ylab="Mean H'",
           xlim=c(-0.1,3),
           data=d2
           )
trellis.device(device="pdf",file="MeanHprimePerLocusLargeEffect.pdf",height=10,width=10)
trellis.par.set("fontsize",list(text=18))
print(p)
dev.off()
