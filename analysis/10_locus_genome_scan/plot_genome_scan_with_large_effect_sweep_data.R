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
    group_by(generation,window,dist,opt,mu) %>%
    summarise_all(funs(mean)) %>% 
    mutate(scaled_time = (generation-50000)/5000) 
soft_only = dbt %>%
    select(-c(repid, locus)) %>%
    filter(nhard==0,nsoft>0,opt>0.1) %>%
    mutate(dist = abs(window-5)) %>%
    group_by(generation,window,dist,opt,mu) %>%
    summarise_all(funs(mean)) %>% 
    mutate(scaled_time = (generation-50000)/5000) 

hard_only=collect(hard_only)
soft_only=collect(soft_only)

COLORS=rev(viridis(length(unique(as.factor(hard_only$dist)))))
KEY=list(space="top",columns=3,title="Distance from window with causal mutations.",
         cex.title=1,points=FALSE,lines=TRUE,just=0.5)
STRIP=strip.custom(strip.names = TRUE,sep=" = ", 
                   var.name = c(expression(mu),expression(z[o])),bg=c("white"))

save_image <- function(stat,img)
{
    trellis.device(device="pdf",file=paste(stat,".pdf",sep=""),height=10,width=10)
    trellis.par.set("fontsize",list(text=18))
    print(img)
    dev.off()
}

X=rbind(hard_only,soft_only)
X=X%>%
    mutate(type=ifelse(nhard>0,'New mutation','Standing var.'))
tajdPlot = xyplot(tajd ~ scaled_time| as.factor(mu)*as.factor(opt):as.factor(type),group=dist,
                  type='l',data=X,
                  par.settings=simpleTheme(col=COLORS),
                  auto.key=KEY,
                  xlab="Time since optimum shift (units of N generations)",
                  ylab="Mean Tajima's D",
                  scales=list(cex=1,alternating=F),
                  strip=STRIP
                  )
save_image('MeanTajdTenLociLargeEffectOnly',tajdPlot)

hprimePlot = xyplot(hprime ~ scaled_time| as.factor(mu)*as.factor(opt):as.factor(type),group=dist,
                  type='l',data=X,
                  par.settings=simpleTheme(col=COLORS),
                  auto.key=KEY,
                  xlab="Time since optimum shift (units of N generations)",
                  ylab="Mean H'",
                  scales=list(cex=1,alternating=F),
                  strip=STRIP)
save_image('MeanHprimeTenLociLargeEffectOnly',hprimePlot)

thetapiPlot = xyplot(thetapi ~ scaled_time| as.factor(mu)*as.factor(opt):as.factor(type),group=dist,
                  type='l',data=X,
                  par.settings=simpleTheme(col=COLORS),
                  auto.key=KEY,
                  xlab="Time since optimum shift (units of N generations)",
                  ylab=expression(paste("Mean ",hat(theta)[pi])),
                  scales=list(cex=1,alternating=F),
                  strip=STRIP)
save_image('MeanPiTenLociLargeEffectOnly',thetapiPlot)

H1Plot = xyplot(H1 ~ scaled_time| as.factor(mu)*as.factor(opt):as.factor(type),group=dist,
                  type='l',data=X,
                  par.settings=simpleTheme(col=COLORS),
                  auto.key=KEY,
                  xlab="Time since optimum shift (units of N generations)",
                  ylab=expression(paste("Mean ",H[1])),
                  scales=list(cex=1,alternating=F),
                  strip=STRIP,xlim=c(-0.05,0.1))
save_image('MeanH1TenLociLargeEffectOnly',H1Plot)

H12Plot = xyplot(H12 ~ scaled_time| as.factor(mu)*as.factor(opt):as.factor(type),group=dist,
                  type='l',data=X,
                  par.settings=simpleTheme(col=COLORS),
                  auto.key=KEY,
                  xlab="Time since optimum shift (units of N generations)",
                  ylab=expression(paste("Mean ",H[12])),
                  scales=list(cex=1,alternating=F),
                  strip=STRIP,xlim=c(-0.05,0.1))
save_image('MeanH12TenLociLargeEffectOnly',H12Plot)

H2H1Plot = xyplot(H2H1 ~ scaled_time| as.factor(mu)*as.factor(opt):as.factor(type),group=dist,
                  type='l',data=X,
                  par.settings=simpleTheme(col=COLORS),
                  auto.key=KEY,
                  xlab="Time since optimum shift (units of N generations)",
                  ylab=expression(paste("Mean ",H[12])),
                  scales=list(cex=1,alternating=F),
                  strip=STRIP,xlim=c(-0.05,0.15))
save_image('MeanH2H1TenLociLargeEffectOnly',H2H1Plot)

mean_nSLPlot = xyplot(mean_nSL ~ scaled_time| as.factor(mu)*as.factor(opt):as.factor(type),group=dist,
                  type='l',data=X,
                  par.settings=simpleTheme(col=COLORS),
                  auto.key=KEY,
                  xlab="Time since optimum shift (units of N generations)",
                  ylab=expression(paste("Mean ",n[SL])),
                  scales=list(cex=1,alternating=F),
                  strip=STRIP,xlim=c(-2,4))
save_image('MeannSLTenLociLargeEffectOnly',mean_nSLPlot)


