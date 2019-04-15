library(dplyr)
library(viridis)
library(DBI)
library(stringr)
library(lattice)

TWINDOW=1e3/11.0
process_files <- function(d, param)
{
    df=NA
    for (i in d)
    {
    print(i)
    x = str_split(i,regex('\\.'))[[1]]
    mu = as.numeric(str_replace(paste(x[2],x[3],sep="."),'mu',''))
    plarge = as.numeric(str_replace(paste(x[4],x[5],sep="."),'plarge',''))
    ghat = 2*sqrt(2)*sqrt(mu)
    db = DBI::dbConnect(RSQLite::SQLite(),i)
    dbt=tbl(db,'data') %>% 
        mutate(dist=abs(window-5)) %>%
        group_by(generation,dist) %>%
        summarise(
                  mD=mean(tajd),
                  mHp=mean(hprime),
                  mz=mean(meanz),
                  mH1=mean(H1),
                  mH12=mean(H12),
                  mH2H1=mean(H2H1))
                  
    results=collect(dbt) %>% mutate(mu=mu,plarge=plarge)
    if(is.na(df))
    {
        df=results
    } else
    {
        df=rbind(df,results)
    }
    DBI::dbDisconnect(db)
}    
    df=df%>% mutate(scaled_time=(generation-5e4),
                    des=param)
    df
}
gamma_small_files = dir('.',pattern='gscan_gamma_small\\..+sqlite3')
gamma_files = dir('.',pattern='gscan_gamma\\..+sqlite3')
gauss_files = dir('.',pattern='gscan\\..+sqlite3')

gauss_df = process_files(gauss_files, 'gauss')
gamma_df = process_files(gamma_files, 'gamma1')
gamma2_df = process_files(gamma_small_files, 'gamma2')

df = rbind(gauss_df, gamma_df, gamma2_df) 

d0 = df %>% filter(dist == 0)

# Plot stats for the middle window
COLORS=rev(viridis(3,alpha=0.5))
STRIP=strip.custom(strip.names = TRUE,sep=" = ", 
                   var.name = c(expression(mu),expression(paste("Pr(|",gamma,"| >= ",hat(gamma),")"))),bg=c("white"))
key_text=list(c(expression(paste(Gamma,", shape = 1.0")),
              expression(paste(Gamma,", shape = 0.5")),
              expression(paste("gaussian"))))
KEY=list(space="top",columns=3,
         title="Distribution of effect sizes",
         cex.title=1,#points=FALSE,
         lines=list(lwd=rep(3,length(COLORS)),col=COLORS),
         just=0.5,
         text=key_text)

p = xyplot(mD~scaled_time|as.factor(mu)*as.factor(plarge),
           data=d0,
           type='l',
           par.settings=simpleTheme(col=COLORS),
           scales=list(cex=1,alternating=F,x=list(rot=45)),
           strip=STRIP,
           group=des,
           lwd=3,
           xlim=c(-0.1,4)*5e3,
           xlab="Generations since optimum shift",
           ylab="Mean Tajima's D in central window",
           key=KEY)
trellis.device('pdf',file='MeanTajimasDVaryPlargeDist0.pdf',height=10,width=10)

trellis.par.set("fontsize",list(text=18))
print(p)
dev.off()

p = xyplot(mHp~scaled_time|as.factor(mu)*as.factor(plarge),
           data=d0,
           type='l',
           par.settings=simpleTheme(col=COLORS),
           scales=list(cex=1,alternating=F,x=list(rot=45)),
           strip=STRIP,
           group=des,
           lwd=3,
           xlim=c(-0.1,4)*5e3,
           xlab="Generations since optimum shift",
           ylab="Mean H' in central window",
           key=KEY)
trellis.device('pdf',file='MeanHprimeVaryPlargeDist0.pdf',height=10,width=10)

trellis.par.set("fontsize",list(text=18))
print(p)
dev.off()

#q("no")

UDIST=unique(as.factor(df$dist))
ICOLORS=viridis(length(UDIST))
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
         text=list(as.character(sort(UDIST))))

# the nSL plots skip distance == 5, so 
# we'll make new colors/key for them:
NSLCOLORS = COLORS[2:6]
NSLKEY=list(space="top",columns=3,
         title="Distance from window with causal mutations.",
         cex.title=1,#points=FALSE,
         lines=list(lwd=rep(3,length(NSLCOLORS)),col=rev(NSLCOLORS)),
         just=0.5,
         text=list(as.character(sort(UDIST)[1:5])))
STRIP=strip.custom(strip.names = TRUE,sep=" = ", 
                   var.name = c(expression(mu),expression(paste("Pr(|",gamma,"| >= ",hat(gamma),")"))),bg=c("white"))
save_image <- function(stat,img)
{
    trellis.device(device="pdf",file=paste(stat,".pdf",sep=""),height=10,width=10)
    trellis.par.set("fontsize",list(text=18))
    print(img)
    dev.off()
}

plotD <- function(df, name)
{
    Dplot = xyplot(mD~scaled_time|as.factor(mu)*as.factor(plarge),data=df,
                   group=factor(dist,levels=rev(sort(unique(df$dist)))),
                   type='l',
                   par.settings=simpleTheme(col=COLORS),
                   key=KEY,
                   lwd=3,
                   scales=list(cex=1,alternating=F,x=list(rot=45)),
                   xlab="Generations since optimum shift",
                   ylab="Mean Tajima's D",
                   strip=STRIP)

    save_image(name,Dplot)
}
plotD(gamma_df,"MeanTajimasDVaryPlargeGamma")
plotD(gamma2_df,"MeanTajimasDVaryPlargeGammaSmall")
plotD(gauss_df,"MeanTajimasDVaryPlargeGaussian")


plotH <- function(df, name)
{
    Hplot = xyplot(mHp~scaled_time|as.factor(mu)*as.factor(plarge),data=df,
                   group=factor(dist,levels=rev(sort(unique(df$dist)))),
                   type='l',
                   key=KEY,
                   par.settings=simpleTheme(col=COLORS),
                   lwd=3,
                   scales=list(cex=1,alternating=F,x=list(rot=45)),
                   xlab="Generations since optimum shift",
                   ylab="Mean H'",
                   strip=STRIP)
    save_image(name,Hplot)
}
plotH(gamma_df,"MeanHprimeVaryPlargeGamma")
plotH(gamma2_df,"MeanHprimeVaryPlargeGammaSmall")
plotH(gauss_df,"MeanHprimeVaryPlargeGaussian")

q("no")

    H1plot = xyplot(mH1~scaled_time|as.factor(mu)*as.factor(plarge),data=df,
               group=factor(dist,levels=rev(sort(unique(df$dist)))),
               type='l',
               par.settings=simpleTheme(col=COLORS),
               key=KEY,
               lwd=3,
               scales=list(cex=1,alternating=F,x=list(rot=45)),
               xlab="Generations since optimum shift",
               ylab=expression(paste("Mean ",H[1])),
               xlim=c(-0.025,0.1)*5e3,
               strip=STRIP)
save_image("H1VaryPlargeGamma",H1plot)
H12plot = xyplot(mH12~scaled_time|as.factor(mu)*as.factor(plarge),data=df,
               group=factor(dist,levels=rev(sort(unique(df$dist)))),
               type='l',
               par.settings=simpleTheme(col=COLORS),
               key=KEY,
               lwd=3,
               scales=list(cex=1,alternating=F,x=list(rot=45)),
               xlab="Generations since optimum shift",
               ylab=expression(paste("Mean ",H[12])),
               xlim=c(-0.025,0.1)*5e3,
               strip=STRIP)
save_image("H12VaryPlargeGamma",H12plot)
H2H1plot = xyplot(mH2H1~scaled_time|as.factor(mu)*as.factor(plarge),data=df,
               group=factor(dist,levels=rev(sort(unique(df$dist)))),
               type='l',
               par.settings=simpleTheme(col=COLORS),
               key=KEY,
               lwd=3,
               scales=list(cex=1,alternating=F,x=list(rot=45)),
               xlab="Generations since optimum shift",
               ylab=expression(paste("Mean ",H[2],'/',H[1])),
               xlim=c(-0.025,0.1)*5e3,
               strip=STRIP)
save_image("H2H1VaryPlargeGamma",H2H1plot)
Zplot = xyplot(mz~scaled_time|as.factor(mu)*as.factor(plarge),data=subset(df,dist<5),
               group=factor(dist,levels=rev(sort(unique(subset(df,dist<5)$dist)))),
               type='l',
               par.settings=simpleTheme(col=NSLCOLORS),
               key=NSLKEY,
               lwd=3,
               scales=list(cex=1,alternating=F,x=list(rot=45)),
               strip=STRIP,
               xlab="Time since optimum shift (units of N generations)",
               ylab="Mean z-score",
               xlim=c(-1.,2)*5e3)
save_image("ZscoreVaryPlargeGamma",Zplot)
