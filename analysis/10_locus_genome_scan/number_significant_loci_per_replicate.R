# The mean number of loci with
# at least one window significant
# at level alpha, where we vary alpha
library(tidyr)
library(dplyr)
library(tibble)
library(readr)
library(stringr)
library(lattice)
library(viridis)

ndist_db = src_sqlite("nulldist.db",create=F)
ndist_data_table = tbl(ndist_db, 'data')

lower05 <- function(x)
{
    return(as.numeric(quantile(x,0.05)))
}
lower01 <- function(x)
{
    return(as.numeric(quantile(x,0.01)))
}
lower001 <- function(x)
{
    return(as.numeric(quantile(x,0.001)))
}
raw_dist = collect(ndist_data_table)
critical_vals = raw_dist %>% summarise_all(funs(lower05,lower01,lower001))

files <- dir(path="../../mlocus_pickle/",pattern="*.genome_scan.db")

data = data.frame()
for (infile in files)
{
    print(infile)
    dbname = paste("../../mlocus_pickle/",infile,sep="")
    indb = src_sqlite(dbname)
    indb_table = tbl(indb,'data')

    fs <- str_split(infile,regex('\\.'))[[1]]
    fs[2]<-str_replace(fs[2],'mu','')
    fs[4]<-str_replace(fs[4],'opt','')
    mu <- as.numeric(paste(fs[2],fs[3],sep='.'))
    opt <- as.numeric(paste(fs[4],fs[5],sep='.'))

    raw_data = collect(indb_table)
    # Is there at least one window significant at alpha?
    sigwin_query = raw_data %>%        
         group_by(generation,locus,repid) %>%
        summarise(tajd05w = ifelse(sum(tajd <= critical_vals$tajd_lower05)>0,1,0),
               tajd01w = ifelse(sum(tajd <= critical_vals$tajd_lower01)>0,1,0),
               tajd001w = ifelse(sum(tajd <= critical_vals$tajd_lower001)>0,1,0),
               hprime05w = ifelse(sum(hprime <= critical_vals$hprime_lower05)>00,1,0),
               hprime01w = ifelse(sum(hprime <= critical_vals$hprime_lower01)>0,1,0),
               hprime001w = ifelse(sum(hprime <= critical_vals$hprime_lower001)>0,1,0))

    #Get how many loci have at least one sig window per rep
    sigwin = collect(sigwin_query)
    sigwin = sigwin %>%
        group_by(repid,generation) %>%
        summarise(tajd05l = sum(tajd05w),
                  tajd01l = sum(tajd01w),
                  tajd001l = sum(tajd001w),
                  hprime05l = sum(hprime05w),
                  hprime01l = sum(hprime01w),
                  hprime001l = sum(hprime001w)) %>%
        group_by(generation) %>%
        summarise(tajd05 = mean(tajd05l),
                  tajd01 = mean(tajd01l),
                  tajd001 = mean(tajd001l),
                  hprime05 = mean(hprime05l),
                  hprime01 = mean(hprime01l),
                  hprime001 = mean(hprime001l)) %>%
        mutate(mu=mu,opt=opt)
    sigwin = sigwin %>% gather(stat, value, tajd05:hprime001) %>% extract(stat, c('stat', 'alpha'), '([^0-9]+)([0-9]+)') %>% mutate(alpha=as.numeric(paste0('.', alpha)))

    print(sigwin)
    data = bind_rows(data,sigwin)
}

KEY=list(space="top",title="Per-window significance threshold",
         cex.title=1,points=FALSE,lines=TRUE,just=0.5)
COLORS=viridis(length(unique(as.factor(data$alpha))))

print(COLORS)
STRIP=strip.custom(strip.names = TRUE,sep=" = ", 
                   var.name = c(expression(mu),expression(z[o])),bg=c("white"))
data=data%>%mutate(scaled_time=(generation-5e4)/5e3)
p = xyplot(value~scaled_time | as.factor(mu)*as.factor(opt),
           data=subset(data,str_detect(stat,"tajd")==T),type='l',
           group=as.factor(alpha),
          par.settings=simpleTheme(col=COLORS,lwd=3),
           strip=STRIP,
          auto.key=KEY,
          xlim=c(-0.5,4),
          xlab="Time since optimum shift (units of N generations)",
          ylab=expression("Mean number of loci with at least one significant window"),
          scales=list(alternating=F))

trellis.device(device="pdf",file="MeanNoSigLociTajD.pdf",height=10,width=10)
trellis.par.set("fontsize",list(text=18))
print(p)
dev.off()
p = xyplot(value~scaled_time | as.factor(mu)*as.factor(opt),
           data=subset(data,str_detect(stat,"tajd")==F),type='l',
           group=as.factor(alpha),
          par.settings=simpleTheme(col=COLORS,lwd=3),
           strip=STRIP,
          auto.key=KEY,
          xlim=c(-0.5,4),
          xlab="Time since optimum shift (units of N generations)",
          ylab=expression("Mean number of loci with at least one significant window"),
          scales=list(alternating=F))

trellis.device(device="pdf",file="MeanNoSigLociHprime.pdf",height=10,width=10)
trellis.par.set("fontsize",list(text=18))
print(p)
dev.off()
