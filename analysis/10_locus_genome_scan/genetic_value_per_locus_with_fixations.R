library(stringr)
library(readr)
library(tidyr)
library(dplyr)
library(lattice)
library(viridis)

process_gvalues <- function (filename, fixations)
{
    fs <- str_split(filename,regex('\\.'))[[1]]
    fs[2]<-str_replace(fs[2],'mu','')
    fs[4]<-str_replace(fs[4],'opt','')
    dbmu <- as.numeric(paste(fs[2],fs[3],sep='.'))
    dbopt <- as.numeric(paste(fs[4],fs[5],sep='.'))

    f = fixations %>% filter(opt==dbopt,mu==dbmu)
    dbname = paste("../../mlocus_pickle/",filename,sep="")
    if(!file.exists(dbname))
    {
        stop("database not found")
    } 
    db = DBI::dbConnect(RSQLite::SQLite(),dbname)
    dbt = tbl(db,'data')

    g = collect(dbt) %>%
        full_join(f,by=c('repid','locus')) %>%
        # The full join will put NA values
        # for all loci that had no sweeps at all.
        # We need to recode those with 0.
        # The join also mangles mu/opt for loci w/no sweeps,
        # so we need to put that back in.
        replace_na(list(h100=0,h200=0,h=0,soft=0,ttl_sweeps=0,opt = dbopt, mu = dbmu)) %>%
        # replace_na changes our int to numeric
        mutate(h100 = as.integer(h100),h200 = as.integer(h200), h = as.integer(h),
               soft = as.integer(soft), ttl_sweeps = as.integer(ttl_sweeps))
    DBI::dbDisconnect(db)
    g
}

fixations = read_delim("raw_fixations.txt.gz",delim=" ")

#See note in file header above re: why we remove none here
fixations = fixations %>% filter(sweep_type != 'none') %>%
    group_by(repid,locus,mu,opt) %>%
    summarise(h100 = as.integer(length(which(sweep_type == 'hard' & g <= 5e4 + 100))),
              h200 = as.integer(length(which(sweep_type == 'hard' & g > 5e4 + 100 & g <= 5e4 + 200))),
              h = as.integer(length(which(sweep_type == 'hard' & g > 5e4 + 200))),
              soft = as.integer(length(which(sweep_type == 'soft')))) %>%
    mutate(ttl_sweeps = h100 + h200 + h + soft)

files <- dir(path="../../mlocus_pickle/",pattern="*.genetic_values_per_locus.db")

d = data.frame() 
d2 = data.frame()
for (i in files)
{
    x = process_gvalues(i, fixations)
    soft = subset(x, soft > 0  & ttl_sweeps == soft) %>%
        rename(softg=g) %>%
        gather(key=g,value=value,softg)
    h100 = subset(x, h100 == 1 & ttl_sweeps == h100) %>%
        rename(h100g=g) %>%
        gather(key=g,value=value,h100g)
    h200 = subset(x, h200 == 1 & ttl_sweeps == h200) %>%
        rename(h200g=g) %>%
        gather(key=g,value=value,h200g)
    h = subset(x, h == 1 & ttl_sweeps == h) %>%
        rename(hg=g) %>%
        gather(key=g,value=value,hg)
    nosweeps = subset(x, ttl_sweeps == 0) %>%
        rename(nosweepsg = g) %>%
        gather(key=g,value=value,nosweepsg)
    somesweeps = subset(x, ttl_sweeps > 0) %>%
        rename(somesweepsg = g) %>%
        gather(key=g,value=value,somesweepsg)
    d = rbind(d,soft,h100,h200,h,nosweeps)
    d2 = rbind(d2,somesweeps,nosweeps)
}

d = d %>%
    group_by(opt,mu,generation,g) %>%
    summarise(mean_g_per_locus = mean(value)) %>%
    mutate(scaled_time = (generation-5e4)/5e4)

d2 = d2 %>%
    group_by(opt,mu,generation,g) %>%
    summarise(mean_g_per_locus = mean(value)) %>%
    mutate(scaled_time = (generation-5e4)/5e4)


KEY=list(space="top",columns=3,title="",
         cex.title=1,points=FALSE,lines=TRUE,just=0.5)
STRIP=strip.custom(strip.names = TRUE, 
                   var.name = c(expression(z[o]),expression(mu)),bg=c("white"))
COLORS=viridis(length(unique(as.factor(d$g))))
p = xyplot(mean_g_per_locus ~ scaled_time|as.factor(opt)*as.factor(mu),data=d,
           type='l',
           lwd=3,
           group = g, par.settings=simpleTheme(col=COLORS,lwd=3),
          auto.key=KEY, xlab="Time since optimum shift (units of N generations)",
          ylab="Mean value of statistic.",
          scales=list(cex=1,alternating=F),
          strip=STRIP,xlim=c(-0.01,0.05))

trellis.device(device="png",file="test.png")
print(p)
dev.off()

KEY=list(space="top",title="",
         cex.title=1,points=FALSE,lines=TRUE,just=0.5)
COLORS=viridis(length(unique(as.factor(d2$g))))
p = xyplot(mean_g_per_locus ~ scaled_time|as.factor(opt)*as.factor(mu),data=d2,
           type='l',
           lwd=3,
           group = g, par.settings=simpleTheme(col=COLORS,lwd=3),
          auto.key=KEY, xlab="Time since optimum shift (units of N generations)",
          ylab="Mean value of statistic.",
          scales=list(cex=1,alternating=F),
          strip=STRIP,xlim=c(-0.01,0.05))

trellis.device(device="png",file="test2.png")
print(p)
dev.off()
