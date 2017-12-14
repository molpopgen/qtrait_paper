# Notes: sweep_type = none
# means fixation prior to optimum shift

# The current problem in this script is
# that some loci have hard and soft sweeps.
# It is not clear which one is getting assigned
# in the join.
library(readr)
library(dplyr)
library(dbplyr)
library(stringr)
library(viridis)
library(lattice)

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

    dt = collect(dbt) %>%
        full_join(f,by=c('repid','locus'))

    print(head(dt))
    dt$sweep_type[which(is.na(dt$sweep_type))] = "none"
    print(unique(dt$sweep_type))

    h = dt %>% filter(sweep_type == 'hard')
    nh = dt %>% filter(sweep_type != 'hard')

    h$sweep_type[which(h$g <= 5e4+100)] = 'hard_recent'
    h$sweep_type[which(h$g > 5e4+100 & h$g <= 5e4+200)] = 'hard_less_recent'
    
    dt = rbind(h,nh)
    dt
}

fixations = read_delim("raw_fixations.txt.gz",delim=" ")
fixations = fixations %>% filter(sweep_type != 'none')

files <- dir(path="../../mlocus_pickle/",pattern="*.genome_scan.db")

data = data.frame()

x=process_genome_scan(files[3],fixations)
x=x%>%
    group_by(sweep_type,generation) %>%
    summarise(md=mean(tajd))

p=xyplot(md~generation,data=x,group=sweep_type,auto.key=TRUE)
trellis.device(device="pdf",file="test.pdf")
print(p)
dev.off()
# for (i in files)
# {
#     data = bind_rows(data,process_genome_scan(i, fixations))
# }
