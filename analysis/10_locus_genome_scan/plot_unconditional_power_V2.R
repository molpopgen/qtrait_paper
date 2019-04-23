# What we want here is to get things like
# Pr(significant test | sweep from new mutation).
# So, this is just Bayes' theorem, meaning we need to get
# Pr(sweep type | significant) Pr(significant)/Pr(sweep type)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(lattice)
library(viridis)

# Read in our null distribution 
ndist_db=src_sqlite('nulldist.db')
ndist_t=tbl(ndist_db,'data')
ndist_raw=collect(ndist_t)

# Because we are going to turn raw window into distance of a 
# window from the window containing selected sites, the quantiles
# corresponding to critical values differ by a factor of 2.
ndist = tibble("dcrit"=c(as.numeric(quantile(ndist_raw$tajd,0.05/(110)))),
               "dcrit2"=c(as.numeric(quantile(ndist_raw$tajd,0.05/(2*110)))),
               "hcrit"=c(as.numeric(quantile(ndist_raw$hprime,0.05/(110)))),
               "hcrit2"=c(as.numeric(quantile(ndist_raw$hprime,0.05/(2*110)))))

# OK, now we have to do heavy lifting, and process the raw genome
# scan data
files <- dir(path="../../mlocus_pickle/",pattern="*.genome_scan.db")
genome_scan_data = NA
for (f in files)
{
    print(f)
    fs <- str_split(f,regex('\\.'))[[1]]
    fs[2]<-str_replace(fs[2],'mu','')
    fs[4]<-str_replace(fs[4],'opt','')
    dbmu <- as.numeric(paste(fs[2],fs[3],sep='.'))
    dbopt <- as.numeric(paste(fs[4],fs[5],sep='.'))
    dbname = paste("../../mlocus_pickle/",f,sep="")
    db <- src_sqlite(dbname)
    dbt <- tbl(db,'data')
    # Here, we add a 0/1 if D and H are more 
    # negative than expected under the null.
    # The multiple test correction is applied
    # here, and it is based on 'dist'.
    dt = collect(dbt) %>% mutate(opt=dbopt,mu=dbmu,dist=abs(window-5),
                                 adjustment=ifelse(dist==0,1,2),
                                 Dsig = ifelse(adjustment == 1,ifelse(tajd<=ndist$dcrit,1,0),ifelse(tajd<=ndist$dcrit2,1,0)),
                                 Hsig = ifelse(adjustment == 1,ifelse(hprime<=ndist$hcrit,1,0),ifelse(hprime<=ndist$hcrit2,1,0))) %>%
    select(repid,locus,window,generation,dist,adjustment,opt,mu,Dsig,Hsig)
    if(is.na(genome_scan_data))
    {
        genome_scan_data=dt
    }
    else
    {
        genome_scan_data=rbind(genome_scan_data,dt)
    }
}
# these are the unconditional 
# probabilities that a test
# is significant, after multiple 
# testing correction
psig = genome_scan_data %>%
    group_by(mu,opt,generation,dist) %>%
    summarise(pDsig = sum(Dsig)/n(),pHsig = sum(Hsig)/n() ,n=n()) %>%
    mutate(scaled_time=(generation-5e4))

# Set up our plots
save_image <- function(stat,img)
{
    trellis.device(device="pdf",file=paste(stat,".pdf",sep=""),height=10,width=10)
    trellis.par.set("fontsize",list(text=18))
    print(img)
    dev.off()
}

# KEY=list(space="top",columns=3,title="Distance from window with causal mutations.",
#          cex.title=1,points=FALSE,lines=TRUE,just=0.5)
STRIP=strip.custom(strip.names = TRUE,sep=" = ", 
                   var.name = c(expression(mu),expression(z[o])),bg=c("white"))

XLIM=c(-0.5,4)*5e3

ICOLORS=viridis(length(unique(as.factor(psig$dist))))
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
nlines = length(unique(psig$dist))
KEY=list(space="top",columns=3,
         title="Distance from window with causal mutations.",
         cex.title=1,#points=FALSE,
         lines=list(lwd=rep(3,length(COLORS)),col=rev(COLORS)),
         just=0.5,
         text=list(as.character(sort(unique(psig$dist)))))

PLOT = xyplot(pDsig~scaled_time | as.factor(mu)*as.factor(opt),
              lwd=3,data=psig,type='l',
              group=factor(dist,levels=rev(sort(unique(psig$dist)))),
              #group=dist,data=s,type='l',
              xlim=XLIM,
              par.settings=simpleTheme(col=COLORS),
              key=KEY,
              scales=list(cex=1,alternating=F),
              strip=STRIP,
              xlab="Generations since optimum shift",
              ylab=expression(paste("Power at ",alpha," = 0.05")))
save_image("PowerTAJD",PLOT)
PLOT = xyplot(pHsig~scaled_time | as.factor(mu)*as.factor(opt),
              lwd=3,data=psig,type='l',
              group=factor(dist,levels=rev(sort(unique(psig$dist)))),
              #group=dist,data=s,type='l',
              xlim=XLIM,
              par.settings=simpleTheme(col=COLORS),
              key=KEY,
              scales=list(cex=1,alternating=F),
              strip=STRIP,
              xlab="Generations since optimum shift",
              ylab=expression(paste("Power at ",alpha," = 0.05")))
save_image("PowerHPRIME",PLOT)
# p = xyplot(pDsig ~ generation|as.factor(mu)*as.factor(opt), data=psig, group=dist,auto.key=TRUE,type='l')
# trellis.device(device="pdf",file="test.pdf",height=10,width=10)
# trellis.par.set("fontsize",list(text=18))
# print(p)
# dev.off()
