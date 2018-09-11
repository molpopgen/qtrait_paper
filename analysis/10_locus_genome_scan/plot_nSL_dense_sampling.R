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
        # Originate < 100 gens after shift
        filter(g < 5e4 + 100) %>%
        # Fix after shift
        filter(ftime >= 5e4) %>%
        # Is of large effect
        filter(abs(s) >= ghat) %>%
        mutate(sweep_type = ifelse(g < 5e4,"stand","new"),
               locus = as.integer(pos/12)) 
    df = collect(q) %>%
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
    mutate(scaled_time=(generation-5e4)/5e3)

conditional_means = data %>% filter(opt > 0.1 & ((nhard == 1 & nsoft == 0) | (nhard == 0 & nsoft == 1))) %>%
    mutate(sweep_type = ifelse(nsoft == 1,"Standing Var.","New mutation")) %>% 
    group_by(generation, dist, mu, opt, sweep_type) %>%
    summarise(mz = mean(mean_nSLz), mH1 = mean(H1), mH12 = mean(H12), mH2H1 = mean(H2H1)) %>%
    mutate(scaled_time=(generation-5e4)/5e3)

ICOLORS=viridis(length(unique(as.factor(unconditional_means$dist))))
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

save_image <- function(stat,img)
{
    trellis.device(device="pdf",file=paste(stat,".pdf",sep=""),height=10,width=10)
    trellis.par.set("fontsize",list(text=18))
    print(img)
    dev.off()
}

# Unconditional nSL plot

nSL = xyplot(mz ~ scaled_time|as.factor(mu)*as.factor(opt), data = subset(unconditional_means, dist<5),
            group=factor(dist,levels=rev(sort(unique(unconditional_means$dist)))),
            type='l',
            par.settings=simpleTheme(col=COLORS),
            key=KEY, lwd=3,
            xlab="Time since optimum shift (units of N generations)",
            ylab=expression(paste("Mean z-score")),
            scales=list(cex=1,alternating=F),
            xlim=c(-0.5,4),
            strip=STRIP)

save_image("MeannSLTenLoci", nSL)

nSL = xyplot(mz ~ scaled_time|as.factor(mu)*as.factor(opt):as.factor(sweep_type),
            data = subset(conditional_means, dist<5),
            group=factor(dist,levels=rev(sort(unique(conditional_means$dist)))),
            type='l',
            par.settings=simpleTheme(col=COLORS),
            key=KEY, lwd=3,
            xlab="Time since optimum shift (units of N generations)",
            ylab=expression(paste("Mean z-score")),
            scales=list(cex=1,alternating=F),
            xlim=c(-0.5,4),
            strip=STRIP)

save_image("MeannSLTenLociLargeEffectOnly", nSL)
