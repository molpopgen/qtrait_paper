library(dplyr)
library(stringr)
library(lattice)
library(viridis)

# NB: thetapi and hprime have labels swapped in the sqlite3 files!

process_genome_scan <- function(filename)
{
    fs <- str_split(filename,regex('\\.'))[[1]]
    fs[2]<-str_replace(fs[2],'mu','')
    fs[4]<-str_replace(fs[4],'opt','')
    mu <- as.numeric(paste(fs[2],fs[3],sep='.'))
    opt <- as.numeric(paste(fs[4],fs[5],sep='.'))

    dbname = paste("../../mlocus_pickle/",filename,sep="")
    db <- src_sqlite(dbname)
    dbt <- tbl(db,'data')

    q <- dbt %>% 
        select(-c(repid, locus)) %>%
        mutate(dist=abs(window-5)) %>% 
        group_by(generation,window,dist) %>%
        summarise_all(funs(mean)) %>% 
        mutate(mu=mu,opt=opt,scaled_time = (generation-50000))#/5000) 

    dt = collect(q)
    dt
}

save_image <- function(stat,img)
{
    trellis.device(device="pdf",file=paste(stat,".pdf",sep=""),height=10,width=10)
    trellis.par.set("fontsize",list(text=18))
    print(img)
    dev.off()
}


files <- dir(path="../../mlocus_pickle/",pattern="*.genome_scan.db")

data = data.frame()

for (i in files)
{
    data = bind_rows(data,process_genome_scan(i))
}

# COLORS=rev(viridis(length(unique(as.factor(data$dist)))))
# KEY=list(space="top",columns=3,title="Distance from window with causal mutations.",
#          cex.title=1,points=FALSE,lines=TRUE,just=0.5)
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
                   var.name = c(expression(mu),expression(z[o])),bg=c("white"))

tajdPlot = xyplot(tajd ~ scaled_time| as.factor(mu)*as.factor(opt),
                  group=factor(dist,levels=rev(sort(unique(data$dist)))),
                  type='l',data=data,
                  par.settings=simpleTheme(col=COLORS),
                  key=KEY, lwd=3,
                  xlab="Generations since optimum shift",
                  ylab="Mean Tajima's D",
                  xlim=c(-2500,4*5e3),
                  scales=list(cex=0.75,alternating=F,x=list(rot=45)),
                  strip=STRIP)
save_image("MeanTajdTenLoci",tajdPlot)

hprimePlot = xyplot(hprime ~ scaled_time| as.factor(mu)*as.factor(opt),
                  group=factor(dist,levels=rev(sort(unique(data$dist)))),
                  type='l',data=data,
                  par.settings=simpleTheme(col=COLORS),
                  key=KEY, lwd=3,
                  xlab="Generations since optimum shift",
                  ylab="Mean H'",
                  scales=list(cex=0.75,alternating=F,x=list(rot=45)),
                  #scales=list(cex=1,alternating=F),
                  xlim=c(-0.5,4)*5e3,
                  strip=STRIP)
save_image("MeanHprimeTenLoci",hprimePlot)

thetapiPlot = xyplot(thetapi ~ scaled_time| as.factor(mu)*as.factor(opt),
                  group=factor(dist,levels=rev(sort(unique(data$dist)))),
                  type='l',data=data,
                  par.settings=simpleTheme(col=COLORS),
                  key=KEY, lwd=3,
                  xlab="Generations since optimum shift",
                  ylab=expression(paste("Mean ",hat(theta)[pi])),
                  scales=list(cex=0.75,alternating=F,x=list(rot=45)),
                  #scales=list(cex=1,alternating=F),
                  strip=STRIP)
save_image('MeanPiTenLoci',thetapiPlot)

reduction = xyplot(thetapi/(1000/11) ~ scaled_time| as.factor(mu)*as.factor(opt),
                  group=factor(dist,levels=rev(sort(unique(data$dist)))),
                  type='l',data=data,
                  par.settings=simpleTheme(col=COLORS),
                  key=KEY, lwd=3,
                  xlab="Generations since optimum shift",
                  ylab=expression(paste("Mean ",hat(theta)[pi],"/",theta)),
                  scales=list(cex=0.75,alternating=F,x=list(rot=45)),
                  #scales=list(cex=1,alternating=F),
                  xlim=c(-0.5,4)*5e3,
                  strip=STRIP)
save_image('MeanReductionDiversityTenLoci',reduction)

H1Plot = xyplot(H1 ~ scaled_time| as.factor(mu)*as.factor(opt),
                  group=factor(dist,levels=rev(sort(unique(data$dist)))),
                  type='l',data=data,
                  par.settings=simpleTheme(col=COLORS),
                  key=KEY, lwd=3,
                  xlab="Generations since optimum shift",
                  ylab=expression(paste("Mean ",H[1])),
                  scales=list(cex=0.75,alternating=F,x=list(rot=45)),
                  # scales=list(cex=1,alternating=F),
                  strip=STRIP,xlim=c(-0.05,0.1)*5e3)
save_image('MeanH1TenLoci',H1Plot)

H12Plot = xyplot(H12 ~ scaled_time| as.factor(mu)*as.factor(opt),
                  group=factor(dist,levels=rev(sort(unique(data$dist)))),
                  type='l',data=data,
                  par.settings=simpleTheme(col=COLORS),
                  key=KEY, lwd=3,
                  xlab="Generations since optimum shift",
                  ylab=expression(paste("Mean ",H[12])),
                  scales=list(cex=0.75,alternating=F,x=list(rot=45)),
                  # scales=list(cex=1,alternating=F),
                  strip=STRIP,xlim=c(-0.05,0.1)*5e3)
save_image('MeanH12TenLoci',H12Plot)

H2H1Plot = xyplot(H2H1 ~ scaled_time| as.factor(mu)*as.factor(opt),
                  group=factor(dist,levels=rev(sort(unique(data$dist)))),
                  type='l',data=data,
                  par.settings=simpleTheme(col=COLORS),
                  key=KEY, lwd=3,
                  xlab="Generations since optimum shift",
                  ylab=expression(paste("Mean ",H[12])),
                  scales=list(cex=0.75,alternating=F,x=list(rot=45)),
                  #scales=list(cex=1,alternating=F),
                  strip=STRIP,xlim=c(-0.05,0.15)*5e3)
save_image('MeanH2H1TenLoci',H2H1Plot)

mean_nSLPlot = xyplot(mean_nSL ~ scaled_time| as.factor(mu)*as.factor(opt),
                  group=factor(dist,levels=rev(sort(unique(data$dist)))),
                  type='l',data=data,
                  par.settings=simpleTheme(col=COLORS),
                  key=KEY, lwd=3,
                  xlab="Generations since optimum shift",
                  ylab=expression(paste("Mean ",n[SL])),
                  scales=list(cex=0.75,alternating=F,x=list(rot=45)),
                  #scales=list(cex=1,alternating=F),
                  xlim=c(-0.5,4)*5e3,
                  strip=STRIP)
save_image('MeannSLTenLoci',mean_nSLPlot)

# This one is garbage:
max_abs_nSLPlot = xyplot(max_abx_nSL ~ scaled_time| as.factor(mu)*as.factor(opt),
                  group=factor(dist,levels=rev(sort(unique(data$dist)))),
                  type='l',data=data,
                  par.settings=simpleTheme(col=COLORS),
                  key=KEY, lwd=3,
                  xlab="Generations since optimum shift",
                  ylab=expression(paste("Mean ",H[12])),
                  scales=list(cex=0.75,alternating=F,x=list(rot=45)),
                  #scales=list(cex=1,alternating=F),
                  strip=STRIP,xlim=c(-2,4)*5e3)
