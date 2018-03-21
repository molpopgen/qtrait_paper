library(stringr)
library(readr)
library(tidyr)
library(dplyr)
library(lattice)
library(viridis)

num_unique_loci <- function(x)
{
    print(length(x))
    return(length(unique(x)))
}

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
        filter(generation>=5e4) %>%
        group_by(repid,generation) %>%
        mutate(rank = row_number(desc(abs(g)))) %>%
        #mutate(rank = dense_rank(desc(abs(g)),ties='first')) %>%
        group_by(rank,repid) %>%
        arrange(generation,locus) %>%
        mutate(nunique_loci = as.integer(cummax(as.numeric(factor(locus,levels=unique(locus)))))) %>%
        arrange(generation,locus) %>%
        group_by(generation,rank) %>%
        summarise(mean_num_unique_loci = mean(nunique_loci)) %>%
        mutate(mu=dbmu,opt=dbopt)


    # ranked=subset(g,repid==0) %>%
    # group_by(rank) %>% 
    # arrange(generation,locus) %>%
    # mutate(nunique_loci = as.integer(cummax(as.numeric(factor(locus,levels=unique(locus)))))) %>%
    # arrange(generation,locus)

    # ranked=NA
    # for(r in unique(g$rank))
    # {
    #     for(repid in unique(g$repid))
    #     {
    #         ri=subset(g,rank==r & repid==repid)
    #         # print(nrow(ri))
    #         cr=array()
    #         for (i in 1:nrow(ri))
    #         {
    #             if(i==1){cr[i]=1}
    #             else
    #             {
    #                 # print(i)
    #                 cr[i] = length(unique(ri$locus[1:(i-1)]))
    #             }
    #         }
    #         ri$cr=as.integer(cr)
    #         if(is.na(ranked)){ranked=ri}
    #         else
    #         {
    #             ranked=rbind(ranked,ri)
    #         }
    #     }
    # }
    # ranked = ranked %>% mutate(mu=dbmu,opt=dbopt)
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
# print(warnings())
d = d %>%
    mutate(scaled_time = (generation-5e4)/5e3)

KEY=list(space="top",columns=3,title="",
         cex.title=1,points=FALSE,lines=TRUE,just=0.5)
STRIP=strip.custom(strip.names = TRUE, sep=" = ",
                   var.name = c(expression(mu),expression(z[o])),bg=c("white"))
COLORS=rev(viridis(length(unique(as.factor(d$rank)))))
p = xyplot(mean_num_unique_loci~scaled_time|as.factor(mu)*as.factor(opt),
           type='l',
           group=rank,
           strip=STRIP,
           #auto.key=KEY,
           par.settings=simpleTheme(col=COLORS,lwd=3),
           scales=list(cex=1,alternating=F),
           xlab="Time since optimum shift (units of N generations)",
           ylab="Mean number of distinct loci at rank",
           xlim=c(-0.01,0.2),
           ylim=c(0.0,8),
           data=d)
trellis.device(device="pdf",file="GeneticValuePerLocusRankChanges.pdf",height=10,width=10)
trellis.par.set("fontsize",list(text=18))
print(p)
dev.off()

