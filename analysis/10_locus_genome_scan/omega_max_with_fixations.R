library(dbplyr)
library(dplyr)
library(tidyr)
library(stringr)
library(lattice)
library(viridis)

dbf = src_sqlite("genome_scan_with_large_effect_sweep_counts.sqlite3")
fdbt = tbl(dbf,'data')

# hard_only = dbt %>%
#     select(-c(repid, locus)) %>%
#     filter(nhard>0,nsoft==0,opt>0.1) %>%
#     mutate(dist = abs(window-5)) %>%
#     group_by(generation,window,dist,opt,mu) %>%
#     summarise_all(funs(mean)) %>% 
#     mutate(scaled_time = (generation-50000)/5000) 
# soft_only = dbt %>%
#     select(-c(repid, locus)) %>%
#     filter(nhard==0,nsoft>0,opt>0.1) %>%
#     mutate(dist = abs(window-5)) %>%
#     group_by(generation,window,dist,opt,mu) %>%
#     summarise_all(funs(mean)) %>% 
#     mutate(scaled_time = (generation-50000)/5000) 
# 
# hard_only=collect(hard_only)
# soft_only=collect(soft_only)
# 
# X=rbind(hard_only,soft_only)
# X=X%>%
#     mutate(type=ifelse(nhard>0,'New mutation','Standing var.'))

# Get the omega max numbers
files <- dir(path="../../mlocus_pickle/",pattern="*.omega_max.db")
omax=NA
for (f in files)
{
    print(f)
    fs <- str_split(f,regex('\\.'))[[1]]
    fs[2]<-str_replace(fs[2],'mu','')
    fs[4]<-str_replace(fs[4],'opt','')
    dbmu <- as.numeric(paste(fs[2],fs[3],sep='.'))
    dbopt <- as.numeric(paste(fs[4],fs[5],sep='.'))
    dbname = paste("../../mlocus_pickle/",f,sep="")
    print(dbname)
    print(file.exists(dbname))
    db <- src_sqlite(dbname)
    dbt <- tbl(db,'data')
    dt = collect(dbt) %>% mutate(opt=dbopt,mu=dbmu)
    if(is.na(omax))
    {
        omax=dt
    }
    else
    {
        omax=rbind(omax,dt)
    }
}

#print(head(omax))

fixations=collect(fdbt)
DATA = left_join(fixations,omax,by=c("repid","generation","locus","window","mu","opt"))
print(head(DATA))
hard_only = DATA %>%
    select(-c(repid, locus)) %>%
    filter(nhard>0,nsoft==0,opt>0.1) %>%
    mutate(dist = abs(window-5)) %>%
    group_by(generation,dist,opt,mu) %>%
    summarise_all(funs(mean)) %>% 
    mutate(scaled_time = (generation-50000)/5000) 
soft_only = DATA %>%
    select(-c(repid, locus)) %>%
    filter(nhard==0,nsoft>0,opt>0.1) %>%
    mutate(dist = abs(window-5)) %>%
    group_by(generation,dist,opt,mu) %>%
    summarise_all(funs(mean)) %>% 
    mutate(scaled_time = (generation-50000)/5000) 

print(head(hard_only))

X=rbind(hard_only,soft_only)
X=X%>%
    mutate(type=ifelse(nhard>0,'New mutation','Standing var.'))

# COLORS=rev(viridis(length(unique(as.factor(hard_only$dist)))))
# KEY=list(space="top",columns=3,title="Distance from window with causal mutations.",
#          cex.title=1,points=FALSE,lines=TRUE,just=0.5)
ICOLORS=viridis(length(unique(as.factor(X$dist))))
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
         text=list(as.character(sort(unique(X$dist)))))
STRIP=strip.custom(strip.names = TRUE,sep=" = ", 
                   var.name = c(expression(mu),expression(z[o])),bg=c("white"))

save_image <- function(stat,img)
{
    trellis.device(device="pdf",file=paste(stat,".pdf",sep=""),height=10,width=10)
    trellis.par.set("fontsize",list(text=18))
    print(img)
    dev.off()
}

tajdPlot = xyplot(omegamax ~ scaled_time| as.factor(mu)*as.factor(opt):as.factor(type),
                  group=factor(dist,levels=rev(sort(unique(X$dist)))),
                  type='l',data=X,
                  par.settings=simpleTheme(col=COLORS),
                  key=KEY, lwd=3,
                  xlab="Time since optimum shift (units of N generations)",
                  ylab=expression(paste("Mean ",omega[max])),
                  scales=list(cex=1,alternating=F),
                  strip=STRIP,
                  xlim=c(-0.05,0.15)
                  )
save_image('MeanOmegaMaxTenLociLargeEffectOnly',tajdPlot)
# 
