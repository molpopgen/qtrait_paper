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
        summarise_each(funs(mean)) %>% 
        mutate(mu=mu,opt=opt,scaled_time = (generation-50000)/5000) %>%
        rename(thetapi=hprime,hprime=thetapi)

    dt = collect(q)
    dt
}

save_image <- function(stat,img)
{
    trellis.device(device="pdf",file=paste(stat,".pdf",sep=""),height=10,width=10)
    print(img)
    dev.off()
}


files <- dir(path="../../mlocus_pickle/",pattern="*.genome_scan.db")

data = data.frame()

for (i in files)
{
    data = rbind(data,process_genome_scan(i))
}

COLORS=viridis(length(unique(as.factor(data$dist))))
KEY=list(space="top",columns=3,title="Distance from window with causal mutations.",
         cex.title=1,points=FALSE,lines=TRUE,just=0.5)
STRIP=strip.custom(strip.names = TRUE, 
                   var.name = c(expression(z[o]),expression(mu)),bg=c("white"))

tajdPlot = xyplot(tajd ~ scaled_time| as.factor(opt)*as.factor(mu),group=dist,
                  type='l',data=data,
                  par.settings=simpleTheme(col=COLORS),
                  auto.key=KEY,
                  xlab="Time since optimum shift (units of N generations)",
                  ylab="Mean Tajima's D",
                  scales=list(cex=1),
                  strip=STRIP)
save_image("10_locus_tajd",tajdPlot)

hprimePlot = xyplot(hprime ~ scaled_time| as.factor(opt)*as.factor(mu),group=dist,
                  type='l',data=data,
                  par.settings=simpleTheme(col=COLORS),
                  auto.key=KEY,
                  xlab="Time since optimum shift (units of N generations)",
                  ylab="Mean H'",
                  scales=list(cex=1),
                  strip=STRIP)
save_image("10_locus_hprime",hprimePlot)

thetapiPlot = xyplot(thetapi ~ scaled_time| as.factor(opt)*as.factor(mu),group=dist,
                  type='l',data=data,
                  par.settings=simpleTheme(col=COLORS),
                  auto.key=KEY,
                  xlab="Time since optimum shift (units of N generations)",
                  ylab=expression(paste("Mean ",hat(theta)[pi])),
                  scales=list(cex=1),
                  strip=STRIP)
save_image('10_locus_thetapi',thetapiPlot)

H1Plot = xyplot(H1 ~ scaled_time| as.factor(opt)*as.factor(mu),group=dist,
                  type='l',data=data,
                  par.settings=simpleTheme(col=COLORS),
                  auto.key=KEY,
                  xlab="Time since optimum shift (units of N generations)",
                  ylab=expression("Mean ",H[1]),
                  scales=list(cex=1),
                  strip=STRIP,xlim=c(-0.05,0.1))
save_image('10_locus_H1',H1Plot)

H12Plot = xyplot(H12 ~ scaled_time| as.factor(opt)*as.factor(mu),group=dist,
                  type='l',data=data,
                  par.settings=simpleTheme(col=COLORS),
                  auto.key=KEY,
                  xlab="Time since optimum shift (units of N generations)",
                  ylab=expression(paste("Mean ",H[12])),
                  scales=list(cex=1),
                  strip=STRIP,xlim=c(-0.05,0.1))
save_image('10_locus_H12',H12Plot)

H2H1Plot = xyplot(H2H1 ~ scaled_time| as.factor(opt)*as.factor(mu),group=dist,
                  type='l',data=data,
                  par.settings=simpleTheme(col=COLORS),
                  auto.key=KEY,
                  xlab="Time since optimum shift (units of N generations)",
                  ylab=expression(paste("Mean ",H[12])),
                  scales=list(cex=1),
                  strip=STRIP,xlim=c(-0.05,0.15))
save_image('10_locus_H2H1',H2H1Plot)

mean_nSLPlot = xyplot(mean_nSL ~ scaled_time| as.factor(opt)*as.factor(mu),group=dist,
                  type='l',data=data,
                  par.settings=simpleTheme(col=COLORS),
                  auto.key=KEY,
                  xlab="Time since optimum shift (units of N generations)",
                  ylab=expression(paste("Mean ",H[12])),
                  scales=list(cex=1),
                  strip=STRIP,xlim=c(-2,4))
save_image('10_locus_mean_nSL',mean_nSLPlot)

# This one is garbage:
max_abs_nSLPlot = xyplot(max_abx_nSL ~ scaled_time| as.factor(opt)*as.factor(mu),group=dist,
                  type='l',data=data,
                  par.settings=simpleTheme(col=COLORS),
                  auto.key=KEY,
                  xlab="Time since optimum shift (units of N generations)",
                  ylab=expression(paste("Mean ",H[12])),
                  scales=list(cex=1),
                  strip=STRIP,xlim=c(-2,4))








