# Simple summaries, not 
# broken down by detailed locus properties
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(RSQLite))


muvals = c("0.00025", "0.001", "0.005")
optvals = c("0.1","0.5","1.0")

data = NA
for (m in muvals)
{
    for(o in optvals)
    {
        fn = paste("../pops.mu",m,".opt",o,".genome_scan.db",sep="")
        db <- src_sqlite(fn)
        dbt <- tbl(db,'data')

        query_win <- dbt %>%
            filter(generation >= 8*5000) %>%
            group_by(generation,window) %>%
            summarise_each(funs(mean),-generation,-locus,-repid) %>%
            mutate(opt=as.numeric(o)) %>%
            mutate(mu=as.numeric(m)) %>%
            mutate(dwindow = abs(window-5.0))

        x = collect(query_win)
        if( is.na(data) )
        {
            data = x
        }
        else
        {
            data = rbind(data,x)
        }
    }
}

gz = gzfile("per_locus_genome_scan_means.gz","wb")
write.table(data,gz,row.names=F,quote=F)
close(gz)
