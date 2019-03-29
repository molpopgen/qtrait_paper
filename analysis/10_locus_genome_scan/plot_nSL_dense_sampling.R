library(DBI)
library(dplyr)
library(tidyr)
library(stringr)
library(lattice)
library(viridis)
library(readr)

fixfiles=dir(path="../../mlocus_pickle/",pattern="*.dense_hapstats_fixations.sqlite3")
statfiles=dir(path="../../mlocus_pickle/",pattern="*.dense_hapstats.sqlite3")


get_fixations = function(fname)
{
    x = str_split(fname,regex('\\.'))[[1]]
    mu = as.numeric(str_replace(paste(x[2],x[3],sep="."),'mu',''))
    opt = as.numeric(str_replace(paste(x[4],x[5],sep="."),'opt',''))
    ghat = 2*sqrt(2)*sqrt(mu)
    db = DBI::dbConnect(RSQLite::SQLite(),paste("../../mlocus_pickle/",fname,sep=""))
    dbt = tbl(db,'data')
    q = dbt %>%
        # This filter is for rapid protoyping this script :)
        # filter(repid < 5) %>%
        # Originate < 100 gens after shift
        filter(g < 5e4 + 100) %>%
        # Fix after shift
        filter(ftime >= 5e4) %>%
        # Is of large effect
        mutate(sweep_type = ifelse(g < 5e4,"stand","new"),
               locus = as.integer(pos/12)) 
    df = collect(q) %>%
        filter(abs(s) >= sqrt(100/5e3)) %>%
        group_by(repid,locus) %>%
        summarise(nhard = as.integer(length(which(sweep_type=="new"))),
                  nsoft = as.integer(length(which(sweep_type=="stand"))))
    DBI::dbDisconnect(db)
    df
}

get_stats = function(fname)
{
    db = DBI::dbConnect(RSQLite::SQLite(),paste("../../mlocus_pickle/",fname,sep=""))
    dbt = tbl(db,'data')
    q = dbt %>% 
        # This filter is for rapid protoyping this script :)
        # filter(repid < 5) %>%
        mutate(dist = abs(window-5))
    df = collect(q)
    DBI::dbDisconnect(db)
    df
}

data = data.frame()
for(i in 1:length(fixfiles))
{
    print(paste(fixfiles[i],statfiles[i]))
    x = str_split(fixfiles[i],regex('\\.'))[[1]]
    mu = as.numeric(str_replace(paste(x[2],x[3],sep="."),'mu',''))
    opt = as.numeric(str_replace(paste(x[4],x[5],sep="."),'opt',''))
    fixations=get_fixations(fixfiles[i])
    stats = get_stats(statfiles[i])
    df = full_join(stats,fixations,by=c("repid","locus")) %>%
    replace_na(list(nhard=0,nsoft=0)) %>%
    mutate(mu=mu,opt=opt)
    data = rbind(data,df)
}

unconditional_means = data %>% group_by(generation, dist, mu, opt) %>%
    summarise(mz = mean(mean_nSLz), mH1 = mean(H1), mH12 = mean(H12), mH2H1 = mean(H2H1)) %>%
    mutate(scaled_time=(generation-5e4))

conditional_means = data %>% filter(opt > 0.1 & ((nhard == 1 & nsoft == 0) | (nhard == 0 & nsoft == 1))) %>%
    mutate(sweep_type = ifelse(nsoft == 1,"Standing Var.","New mutation")) %>% 
    group_by(generation, dist, mu, opt, sweep_type) %>%
    summarise(mz = mean(mean_nSLz), mH1 = mean(H1), mH12 = mean(H12), mH2H1 = mean(H2H1)) %>%
    mutate(scaled_time=(generation-5e4))

UDIST=unique(as.factor(unconditional_means$dist))
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
print(COLORS)
print(NSLCOLORS)
NSLKEY=list(space="top",columns=3,
         title="Distance from window with causal mutations.",
         cex.title=1,#points=FALSE,
         lines=list(lwd=rep(3,length(NSLCOLORS)),col=rev(NSLCOLORS)),
         just=0.5,
         text=list(as.character(sort(UDIST)[1:5])))
STRIP=strip.custom(strip.names = TRUE,sep=" = ", 
                   var.name = c(expression(mu),expression(z[o])),bg=c("white"))

save_image <- function(stat,img)
{
    trellis.device(device="pdf",file=paste(stat,".pdf",sep=""),height=10,width=10)
    trellis.par.set("fontsize",list(text=18))
    print(img)
    dev.off()
}

# Unconditional nSL plot

nSL = xyplot(mz ~ scaled_time|as.factor(mu)*as.factor(opt), data = subset(unconditional_means, dist<5),
            group=factor(dist,levels=rev(sort(unique(UDIST)))[1:5]),
            type='l',
            par.settings=simpleTheme(col=NSLCOLORS),
            key=NSLKEY, lwd=3,
            xlab="Generations since optimum shift",
            ylab=expression(paste("Mean z-score")),
            scales=list(cex=0.75,alternating=F,x=list(rot=45),
                        y=list(at=c(-0.05,0,0.05),
                               labels=c('-0.05','0.0','0.05'))),
            xlim=c(-0.5,4)*5e3,ylim=c(-0.06,0.06),
            strip=STRIP)

save_image("MeannSLTenLoci", nSL)

# Conditional nSL plot

nSL = xyplot(mz ~ scaled_time|as.factor(mu)*as.factor(opt):as.factor(sweep_type),
            data = subset(conditional_means, dist<5),
            group=factor(dist,levels=rev(sort(unique(UDIST)))[1:5]),
            type='l',
            par.settings=simpleTheme(col=NSLCOLORS),
            key=NSLKEY, lwd=3,
            xlab="Generations since optimum shift",
            ylab=expression(paste("Mean z-score")),
            scales=list(cex=0.75,alternating=F,x=list(rot=45)),
            xlim=c(-0.5,4)*5e3, ylim=c(-0.25,0.35),
            strip=STRIP)

save_image("MeannSLTenLociLargeEffectOnly", nSL)

# Unconditional plots for H1, H12, and H2/H1

H1plot = xyplot(mH1 ~ scaled_time|as.factor(mu)*as.factor(opt), data = unconditional_means,
            group=factor(dist,levels=rev(sort(unique(UDIST)))),
            type='l',
            par.settings=simpleTheme(col=COLORS),
            key=KEY, lwd=3,
            xlab="Generations since optimum shift",
            ylab=expression(paste("Mean H1")),
            scales=list(cex=0.75,alternating=F,x=list(rot=45)),
            xlim=c(-0.1,0.2)*5e3,
            strip=STRIP)

save_image("MeanH1TenLoci", H1plot)

H12plot = xyplot(mH12 ~ scaled_time|as.factor(mu)*as.factor(opt), data = unconditional_means,
            group=factor(dist,levels=rev(sort(unique(UDIST)))),
            type='l',
            par.settings=simpleTheme(col=COLORS),
            key=KEY, lwd=3,
            xlab="Generations since optimum shift",
            ylab=expression(paste("Mean H12")),
            scales=list(cex=0.75,alternating=F,x=list(rot=45)),
            xlim=c(-0.1,0.2)*5e3,
            strip=STRIP)

save_image("MeanH12TenLoci", H12plot)

H2H1plot = xyplot(mH2H1 ~ scaled_time|as.factor(mu)*as.factor(opt), data = unconditional_means,
            group=factor(dist,levels=rev(sort(unique(UDIST)))),
            type='l',
            par.settings=simpleTheme(col=COLORS),
            key=KEY, lwd=3,
            xlab="Generations since optimum shift",
            ylab=expression(paste("Mean H2H1")),
            scales=list(cex=0.75,alternating=F,x=list(rot=45)),
            xlim=c(-0.1,0.2)*5e3,
            strip=STRIP)

save_image("MeanH2H1TenLoci", H2H1plot)

# Conditional plots for H1, H12, and H2/H1

H1plot = xyplot(mH1 ~ scaled_time|as.factor(mu)*as.factor(opt):as.factor(sweep_type), data = conditional_means,
            group=factor(dist,levels=rev(sort(unique(UDIST)))),
            type='l',
            par.settings=simpleTheme(col=COLORS),
            key=KEY, lwd=3,
            xlab="Generations since optimum shift",
            ylab=expression(paste("Mean H1")),
            scales=list(cex=0.75,alternating=F,x=list(rot=45)),
            xlim=c(-0.1,0.2)*5e3,
            strip=STRIP)

save_image("MeanH1TenLociLargeEffectOnly", H1plot)

H12plot = xyplot(mH12 ~ scaled_time|as.factor(mu)*as.factor(opt):as.factor(sweep_type), data = conditional_means,
            group=factor(dist,levels=rev(sort(unique(UDIST)))),
            type='l',
            par.settings=simpleTheme(col=COLORS),
            key=KEY, lwd=3,
            xlab="Generations since optimum shift",
            ylab=expression(paste("Mean H12")),
            scales=list(cex=0.75,alternating=F,x=list(rot=45)),
            xlim=c(-0.1,0.2)*5e3,
            strip=STRIP)

save_image("MeanH12TenLociLargeEffectOnly", H12plot)

H2H1plot = xyplot(mH2H1 ~ scaled_time|as.factor(mu)*as.factor(opt):as.factor(sweep_type), data = conditional_means,
            group=factor(dist,levels=rev(sort(unique(UDIST)))),
            type='l',
            par.settings=simpleTheme(col=COLORS),
            key=KEY, lwd=3,
            xlab="Generations since optimum shift",
            ylab=expression(paste("Mean H2H1")),
            scales=list(cex=0.75,alternating=F,x=list(rot=45)),
            xlim=c(-0.1,0.2)*5e3,
            strip=STRIP)

save_image("MeanH2H1TenLociLargeEffectOnly", H2H1plot)
