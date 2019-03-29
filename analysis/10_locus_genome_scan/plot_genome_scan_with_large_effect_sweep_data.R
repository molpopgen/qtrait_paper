library(dbplyr)
library(dplyr)
library(tidyr)
library(lattice)
library(viridis)

dbf = src_sqlite("genome_scan_with_large_effect_sweep_counts.sqlite3")
dbt = tbl(dbf,'data')

hard_only = dbt %>%
    select(-c(repid, locus)) %>%
    filter(nhard>0,nsoft==0,opt>0.1) %>%
    mutate(dist = abs(window-5)) %>%
    group_by(generation,dist,opt,mu) %>%
    #group_by(generation,window,dist,opt,mu) %>%
    summarise_all(funs(mean)) %>% 
    mutate(scaled_time = (generation-50000)) 
soft_only = dbt %>%
    select(-c(repid, locus)) %>%
    filter(nhard==0,nsoft>0,opt>0.1) %>%
    mutate(dist = abs(window-5)) %>%
    group_by(generation,dist,opt,mu) %>%
    summarise_all(funs(mean)) %>% 
    mutate(scaled_time = (generation-50000)) 

hard_only=collect(hard_only)
soft_only=collect(soft_only)

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

tajdPlot = xyplot(tajd ~ scaled_time| as.factor(mu)*as.factor(opt):as.factor(type),
                  group=factor(dist,levels=rev(sort(unique(X$dist)))),
                  type='l',data=X,
                  par.settings=simpleTheme(col=COLORS),
                  key=KEY, lwd=3,
                  xlab="Generations since optimum shift",
                  ylab="Mean Tajima's D",
                  scales=list(cex=1,alternating=F,x=list(rot=45)),
                  xlim=c(-0.5,4)*5e3,
                  strip=STRIP
                  )
save_image('MeanTajdTenLociLargeEffectOnly',tajdPlot)

hprimePlot = xyplot(hprime ~ scaled_time| as.factor(mu)*as.factor(opt):as.factor(type),
                  group=factor(dist,levels=rev(sort(unique(X$dist)))),
                  type='l',data=X,
                  par.settings=simpleTheme(col=COLORS),
                  key=KEY, lwd=3,
                  xlab="Generations since optimum shift",
                  ylab="Mean H'",
                  scales=list(cex=1,alternating=F,x=list(rot=45)),
                  xlim=c(-0.5,4)*5e3,
                  strip=STRIP)
save_image('MeanHprimeTenLociLargeEffectOnly',hprimePlot)

thetapiPlot = xyplot(thetapi ~ scaled_time| as.factor(mu)*as.factor(opt):as.factor(type),
                  group=factor(dist,levels=rev(sort(unique(X$dist)))),
                  type='l',data=X,
                  par.settings=simpleTheme(col=COLORS),
                  key=KEY, lwd=3,
                  xlab="Generations since optimum shift",
                  ylab=expression(paste("Mean ",hat(theta)[pi])),
                  scales=list(cex=1,alternating=F,x=list(rot=45)),
                  xlim=c(-0.5,4)*5e3,
                  strip=STRIP)
save_image('MeanPiTenLociLargeEffectOnly',thetapiPlot)

ReductionPlot = xyplot(thetapi/(1000./11.) ~ scaled_time| as.factor(mu)*as.factor(opt):as.factor(type),
                  group=factor(dist,levels=rev(sort(unique(X$dist)))),
                  type='l',data=X,
                  par.settings=simpleTheme(col=COLORS),
                  key=KEY, lwd=3,
                  xlab="Generations since optimum shift",
                  ylab=expression(paste("Mean ",hat(theta)[pi],"/",theta)),
                  scales=list(cex=1,alternating=F,x=list(rot=45)),
                  xlim=c(-0.5,4)*5e3,
                  strip=STRIP)
save_image('MeanReductionDiversityTenLociLargeEffectOnly',ReductionPlot)
