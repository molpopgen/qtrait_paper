library(dbplyr)
library(dplyr)
library(tidyr)
library(lattice)
library(viridis)
library(readr)
library(stringr)
library(RSQLite)

process_genome_scan <- function(filename, fixations)
{
    fs <- str_split(filename,regex('\\.'))[[1]]
    fs[2]<-str_replace(fs[2],'mu','')
    fs[4]<-str_replace(fs[4],'opt','')
    dbmu <- as.numeric(paste(fs[2],fs[3],sep='.'))
    dbopt <- as.numeric(paste(fs[4],fs[5],sep='.'))

    f = fixations %>% filter(opt==dbopt,mu==dbmu)
    dbname = paste("../../mlocus_pickle/",filename,sep="")
    db <- src_sqlite(dbname)
    dbt <- tbl(db,'data')

    gammahat = 2*sqrt(2)*sqrt(dbmu)
    dt = collect(dbt) %>%
        full_join(f,by=c('repid','locus'))  %>%
        na.omit() 
    dt
}
fixations_raw = read_delim("raw_fixations.txt.gz",delim=" ")

fixations = fixations_raw %>%
    filter(sweep_type != 'none',(g-5e4)/5e3<1e2/5e3,
           abs(s) >= 2*sqrt(2)*sqrt(mu)) %>%
    group_by(repid,locus,mu,opt) %>%
    summarise(nhard=length(which(sweep_type=='hard')),nsoft=length(which(sweep_type=='soft')))

ndist_db=src_sqlite('nulldist.db')
ndist_t=tbl(ndist_db,'data')
 
lower05 <- function(x)
{
    return(as.numeric(quantile(x,0.05)))
}
ndist=collect(ndist_t)%>%summarise_all(lower05)
print(head(ndist))

files <- dir(path="../../mlocus_pickle/",pattern="*.genome_scan.db")
HARD=NA
SOFT=NA
#for (i in files[1:1])
for (i in files)
{
    dfi = process_genome_scan(i,fixations) %>%
        mutate(Dsig=as.integer(ifelse(tajd<ndist$tajd,1,0)),
               Hsig=as.integer(ifelse(hprime<ndist$hprime,1,0)))
    if(is.na(HARD))
    {
        HARD=subset(dfi,nhard>0 & nsoft==0)
    }
    else
    { 
        HARD=rbind(HARD,subset(dfi,nhard>0 & nsoft==0))
    }

    if(is.na(SOFT))
    {
        SOFT=subset(dfi,nhard==0 & nsoft>0)
    }
    else
    { 
        SOFT=rbind(SOFT,subset(dfi,nhard==0 & nsoft>0))
    }
}

HARD=HARD %>% mutate(sweep_type='New mutation')
SOFT=SOFT %>% mutate(sweep_type='Standing var.')
data=rbind(HARD,SOFT) %>%
    group_by(sweep_type,mu,opt,generation,window) %>%
    summarise(n=n(),nDsig=sum(Dsig),nHsig=sum(Hsig))%>%
    mutate(scaled_time=(generation-5e4)/5e3,
           dist=as.integer(abs(window-5))) %>%
    mutate(nDsig_multitest = nDsig/(11*256*ifelse(dist==0,1,2)),
           nHsig_multitest = nHsig/(11*256*ifelse(dist==0,1,2))) %>%
    select(-nDsig,-nHsig,-window)

print(head(data))

unconditional_sig = read_delim('nsig_per_window_tidied.txt.gz',delim=" ")
unconditional_sig_tajd = unconditional_sig %>% filter(stat=='tajd',alpha==0.05) %>%
    group_by(opt,mu,generation,window) %>%
    summarise(nDSig_uncond = sum(value)) %>%
    mutate(dist=as.integer(abs(window-5))) %>%
    mutate(nDSig_uncond_multitest = nDSig_uncond/(11*256*ifelse(dist==0,1,2))) %>%
    select(-window,-nDSig_uncond)

unconditional_sig_hprime = unconditional_sig %>% filter(stat=='hprime',alpha==0.05) %>%
    group_by(opt,mu,generation,window) %>%
    summarise(nHSig_uncond = sum(value)) %>%
    mutate(dist=as.integer(abs(window-5))) %>%
    mutate(nHSig_uncond_multitest = nHSig_uncond/(11*256*ifelse(dist==0,1,2))) %>%
    select(-window,-nHSig_uncond)

print(head(unconditional_sig_tajd))

data = data %>% 
    full_join(unconditional_sig_tajd,by=c("opt","mu","generation","dist")) %>%
    full_join(unconditional_sig_hprime,by=c("opt","mu","generation","dist")) %>%
    mutate(fracD = nDsig_multitest/nDSig_uncond_multitest,
           fracH = nHsig_multitest/nHSig_uncond_multitest) %>%
    select(-n,-nDsig_multitest,-nHsig_multitest,-nDSig_uncond_multitest,-nHSig_uncond_multitest) %>%
    filter(opt > 0.1)

# print(head(subset(data,sweep_type=="New mutation")))
# print(head(subset(data,sweep_type=="Standing var.")))
# q("no")

#ICOLORS=rev(viridis(length(unique(as.factor(data$dist)))))
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
nlines = length(unique(data$dist))
KEY=list(space="top",columns=3,title="Distance from window with causal mutations.",
         cex.title=1,#points=FALSE,
         lines=list(lwd=rep(3,length(COLORS)),col=rev(COLORS)),
         #lines=list(lty=rep(1,nlines),lwd=rep(3,nlines),col=rev(COLORS)),
         just=0.5,
         #col = text color...
         # col=rev(COLORS),
         text=list(as.character(sort(unique(data$dist)))))
STRIP=strip.custom(strip.names = TRUE,sep=" = ", 
                   var.name = c(expression(mu),expression(z[o])),bg=c("white"))
#p=xyplot(nDsig/(n*adjustment*11)~scaled_time|as.factor(mu)*as.factor(opt):as.factor(sweep_type),

XLIM=c(-0.5,4)
p=xyplot(fracD~scaled_time|as.factor(mu)*as.factor(opt):as.factor(sweep_type),
         type='l',
         #group=dist,
         group=factor(dist,levels=rev(sort(unique(data$dist)))),
         lwd=3,
         data=data,
         par.settings=simpleTheme(col=COLORS),
         key=KEY,
         xlim=XLIM,
         xlab="Time since optimum shift (units of N generations)",
         ylab="Fraction of significant windows",
         scales=list(cex=1,alternating=F,y=list(at=c(0,0.25,0.5))),
         strip=STRIP#,
       #   panel = function(x, y) {
       #   panel.xyplot(x, y)
       #   panel.abline(h=0.05)
       # }
        )

trellis.device(device='pdf',file='TajDFractionSigLargeEffect.pdf',height=10,width=10)
trellis.par.set("fontsize",list(text=18))
print(p)
dev.off()

p=xyplot(fracH~scaled_time|as.factor(mu)*as.factor(opt):as.factor(sweep_type),
         type='l',
         #group=dist,
         group=factor(dist,levels=rev(sort(unique(data$dist)))),
         data=data,
         lwd=3,
         par.settings=simpleTheme(col=COLORS),
         key=KEY,
         xlim=XLIM,
         xlab="Time since optimum shift (units of N generations)",
         ylab="Fraction of significant windows",
         scales=list(cex=1,alternating=F,y=list(at=c(0,0.25,0.5,0.75))),
         strip=STRIP#,
       #   panel = function(x, y) {
       #   panel.xyplot(x, y)
       #   panel.abline(h=0.05)
       # }
        )

trellis.device(device='pdf',file='HprimeFractionSigLargeEffect.pdf',height=10,width=10)
trellis.par.set("fontsize",list(text=18))
print(p)
dev.off()

