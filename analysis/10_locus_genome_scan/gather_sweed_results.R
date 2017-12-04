library(dplyr)
library(tibble)
library(readr)
library(stringr)

files <- dir(path="../../mlocus_pickle/",pattern="*.sweed.db")

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

    raw_data = collect(indb_table) %>%
        mutate(mu=mu,opt=opt)

    data = bind_rows(data,raw_data)
}

of = gzfile("raw_sweed.txt.gz","wb")
write_delim(data,of)
close(of)
