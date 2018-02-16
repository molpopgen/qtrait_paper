# The mean number of loci with
# at least one window significant
# at level alpha, where we vary alpha
library(dplyr)
library(tibble)
library(readr)
library(stringr)
library(lattice)

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
    sigwin = collect(sigwin_query,n=50)
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
    print(sigwin)
    data = bind_rows(data,sigwin)
}

p = xyplot(tajd05~generation | as.factor(mu)*as.factor(opt),
           data=data,type='l')

trellis.device(device="pdf",file="test.pdf",height=10,width=10)
trellis.par.set("fontsize",list(text=18))
print(p)
dev.off()
