# What we want here is to get things like
# Pr(significant test | sweep from new mutation).
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(lattice)
library(viridis)

# We need to deal with loci that had zero fixations.
# We do so by generating a table of all loci we simulated.
# When we join this table in below, we get NA for any loci 
# w/o fixations of large effect, and we then mutate NA -> 0.
dummy = tibble(mu = c(rep(2.5e-4,256*3*10),rep(1e-3,256*3*10),rep(5e-3,256*3*10)),
               opt = rep(c(rep(0.1,256*10),rep(0.5,256*10),rep(1.0,256*10)),3),
               locus= as.integer(rep(seq(0,9,1),256*3*3)),
               repid=rep(0:255,each=10,times=9))

# Tally number of sweeps per locus, for each replicate and 
# each parameter combo
fixations = read_delim("raw_fixations.txt.gz",delim=" ") %>%
    # filter on origin time of mutation, removing sweeps
    # fixing prior to opt. shift, and filtering on fixations
    # of large effect
    #filter((g-5e4)/5e3<1e2/5e3,sweep_type != 'none',abs(s)>=2*sqrt(2)*sqrt(mu)) %>%
    filter((g-5e4)/5e3<1e2/5e3,sweep_type != 'none',abs(s)>=sqrt(100/5e3)) %>%
    group_by(mu,opt,repid,locus,sweep_type) %>% tally() %>%
    spread(sweep_type,n) %>% #This turns the sweep type into a column and assigns how many sweeps of each type to the locus
    right_join(dummy,by=c("mu","opt","locus","repid")) %>%
    mutate_all(funs(replace(.,is.na(.),0)))

# This holds the probability of large-effect sweeps
# at an individual locus, 
# from new mutations and from standing variation,
# given that they arose w/in 100 gens of opt. shift
# We do not need these #s for this calculation, but
# it is good to know for the paper:
psweep_at_a_locus = fixations %>%
    group_by(mu,opt) %>%
    summarise(phard = sum(ifelse(hard>0,1,0))/n(),
              psoft=sum(ifelse(soft>0,1,0))/n(),
              pboth=sum(ifelse(hard>0&soft>0,1,0))/n(),
              pnone=sum(ifelse(hard==0&soft==0,1,0))/n())

psweep_per_replicate = fixations %>%
    group_by(mu,opt,repid) %>% 
    summarise(replicate_had_hard_sweeps = ifelse(sum(ifelse(hard>0,1,0))>0,1,0),
              replicate_had_soft_sweeps = ifelse(sum(ifelse(soft>0,1,0))>0,1,0),
              replicate_had_hard_and_soft_sweeps = ifelse(replicate_had_hard_sweeps>0&replicate_had_soft_sweeps>0,1,0),
              replicate_had_no_sweeps = ifelse(replicate_had_hard_sweeps==0&replicate_had_soft_sweeps==0,1,0)) %>%
    summarise(phard = sum(replicate_had_hard_sweeps)/n(),
              psoft = sum(replicate_had_soft_sweeps)/n(),
              pboth = sum(replicate_had_hard_and_soft_sweeps)/n(),
              pnone = sum(replicate_had_no_sweeps)/n())

# print(psweep_at_a_locus)
# print(psweep_per_replicate)
# q('no')

# Read in our null distribution 
ndist_db=src_sqlite('nulldist.db')
ndist_t=tbl(ndist_db,'data')
ndist_raw=collect(ndist_t)

# Because we are going to turn raw window into distance of a 
# window from the window containing selected sites, the quantiles
# corresponding to critical values differ by a factor of 2.
ndist = tibble("dcrit"=c(as.numeric(quantile(ndist_raw$tajd,0.05/110))),
               #"dcrit2"=c(as.numeric(quantile(ndist_raw$tajd,0.05/220))),
               "hcrit"=c(as.numeric(quantile(ndist_raw$hprime,0.05/110))))
               #"hcrit2"=c(as.numeric(quantile(ndist_raw$hprime,0.05/220))))

# OK, now we have to do heavy lifting, and process the raw genome
# scan data
files <- dir(path="../../mlocus_pickle/",pattern="*.genome_scan.db")

genome_scan_data = NA
for (f in files)
{
    print(f)
    fs <- str_split(f,regex('\\.'))[[1]]
    fs[2]<-str_replace(fs[2],'mu','')
    fs[4]<-str_replace(fs[4],'opt','')
    dbmu <- as.numeric(paste(fs[2],fs[3],sep='.'))
    dbopt <- as.numeric(paste(fs[4],fs[5],sep='.'))
    dbname = paste("../../mlocus_pickle/",f,sep="")
    db <- src_sqlite(dbname)
    dbt <- tbl(db,'data')
    # Here, we add a 0/1 if D and H are more 
    # negative than expected under the null.
    # The multiple test correction is applied
    # here, and it is based on 'dist'.
    dt = collect(dbt) %>%
        # We only analyze the window w/selected
        # mutations.  The hope is that this simplifies
        # the plotting down stream
        filter(generation < 5e3*14,window==5) %>% # We won't plot more than 4N generations past shift... 
        select(repid,locus,generation,tajd,hprime) %>% # focus on the classics...
        mutate(opt=dbopt,mu=dbmu,#dist=abs(window-5),
               Dsig = ifelse(tajd <= ndist$dcrit,1,0),
               Hsig = ifelse(hprime <= ndist$dcrit,1,0)) %>%
                                 #adjustment=ifelse(dist==0,1,2),
                                 #Dsig = ifelse(adjustment == 1,ifelse(tajd<=ndist$dcrit,1,0),ifelse(tajd<=ndist$dcrit2,1,0)),
                                 #Hsig = ifelse(adjustment == 1,ifelse(hprime<=ndist$hcrit,1,0),ifelse(hprime<=ndist$hcrit2,1,0))) %>%
    select(repid,locus,generation,opt,mu,Dsig,Hsig)
    #select(repid,locus,window,generation,dist,adjustment,opt,mu,Dsig,Hsig)
    if(is.na(genome_scan_data))
    {
        genome_scan_data=dt
    }
    else
    {
        genome_scan_data=rbind(genome_scan_data,dt)
    }
}


# We now need Pr(sweep type and sig), which means some joining...

print(nrow(genome_scan_data))
print(nrow(fixations))
print(names(genome_scan_data))
psweep_and_sig = genome_scan_data %>%
    left_join(fixations,by=c("opt","mu","repid","locus")) %>%
    #filter((hard>0&soft==0)|(hard==0&soft>0)) %>%  #remove that nasty class of both, which is minority of data
    mutate(# Dsig_and_hard = ifelse(Dsig>0&hard>0,1,0),
           # Dsig_and_soft = ifelse(Dsig>0&soft>0,1,0),
           # Hsig_and_hard = ifelse(Hsig>0&hard>0,1,0),
           # Hsig_and_soft = ifelse(Hsig>0&soft>0,1,0),
           # Dsig_and_sweep = ifelse(Dsig>0,1,0),
           # Hsig_and_sweep = ifelse(Hsig>0,1,0),
           sweep_type = case_when((hard > 0 & soft == 0) ~ 0,
                            (hard == 0 & soft > 0) ~ 1,
                            (hard > 0 & soft > 0) ~ 2,
                            (hard+soft == 0) ~ 3)) %>%
           #type = ifelse( hard > 0 & soft == 0, "New mutation",
           #              ifelse(hard == 0 & soft > 0,"Standing var.", "Both"))) %>%
    left_join(psweep_at_a_locus,by=c("opt","mu")) %>%
    filter(sweep_type != 2) # Do not plot category 2, as there are so few replicates.

print(head(psweep_and_sig))
print(psweep_at_a_locus)
print(psweep_and_sig %>%
      group_by(opt,mu,generation) %>% # ,locus) %>%
      summarise(ptot = sum(phard+psoft+pboth+pnone)/n(),n=n()))
# print(psweep_and_sig)
# print(psweep_and_sig%>%filter(Dsig_and_sweep>0)%>%select(repid,locus,dist,type,opt,mu,phard,psoft))
psweep_and_sig = psweep_and_sig %>%
    group_by(opt,mu,generation,sweep_type) %>%
    summarise(pDsig_given_sweep = sum(Dsig)/n(),# (sum(Dsig_and_sweep)/256)/case_when(sweep_type == "New mutation" ~ phard,
                                  #                                      sweep_type == "Standing var." ~ psoft,
                                  #                                      sweep_type == "Both" ~ pboth,
                                  #                                      sweep_type == "None" ~ pnone)[1],
              pHsig_given_sweep = sum(Hsig)/n()) %>%
    mutate(scaled_time = (generation-5e4))


# p = xyplot(pDsig_given_sweep~scaled_time | as.factor(mu)*as.factor(opt),group=sweep_type,data=psweep_and_sig,type='l',auto.key=TRUE)
# #p = xyplot(pDsig_given_hard_sweep~scaled_time | as.factor(mu)*as.factor(opt),group=dist,data=genome_scan_data,type='l')
# trellis.device(device="pdf",file="test.pdf",height=10,width=10)
# trellis.par.set("fontsize",list(text=18))
# print(p)
ICOLORS=viridis(length(unique(as.factor(psweep_and_sig$sweep_type))))
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
nlines = length(unique(psweep_and_sig$sweep_type))
KEY=list(space="top",columns=3,
         title="Type of sweep",
         cex.title=1,#points=FALSE,
         lines=list(lwd=rep(3,length(COLORS)),col=rev(COLORS)),
         just=0.5,
         text=list(c("New mutation","Standing var.","None")))

STRIP=strip.custom(strip.names = TRUE,sep=" = ", 
                   var.name = c(expression(mu),expression(z[o])),bg=c("white"))
XLIM=c(-0.5,4)*5e3
# Set up our plots
save_image <- function(stat,img)
{
    trellis.device(device="pdf",file=paste(stat,".pdf",sep=""),height=10,width=10)
    trellis.par.set("fontsize",list(text=18))
    print(img)
    dev.off()
}
print(unique(factor(psweep_and_sig$sweep_type,levels=sort(unique(psweep_and_sig$sweep_type)))))
# p = xyplot(pDsig_given_sweep~scaled_time | as.factor(mu)*as.factor(opt),group=sweep_type,data=psweep_and_sig,type='l',auto.key=TRUE)
PLOT = xyplot(pDsig_given_sweep~scaled_time | as.factor(mu)*as.factor(opt),
              lwd=3,data=psweep_and_sig,type='l',
              group=factor(sweep_type,levels=rev(sort(unique(sweep_type)))),
              #group=dist,data=s,type='l',
              xlim=XLIM,
              par.settings=simpleTheme(col=COLORS),
              key=KEY,
              #auto.key=TRUE,
              scales=list(cex=1,alternating=F),
              strip=STRIP,
              xlab="Generations since optimum shift",
              ylab=expression("Pr(significant|sweep type)"))
save_image("TajDFractionSigLargeEffect",PLOT)
PLOT = xyplot(pHsig_given_sweep~scaled_time | as.factor(mu)*as.factor(opt),
              lwd=3,data=psweep_and_sig,type='l',
              group=factor(sweep_type,levels=rev(sort(unique(sweep_type)))),
              #group=dist,data=s,type='l',
              xlim=XLIM,
              par.settings=simpleTheme(col=COLORS),
              key=KEY,
              #auto.key=TRUE,
              scales=list(cex=1,alternating=F),
              strip=STRIP,
              xlab="Generations since optimum shift",
              ylab=expression("Pr(significant|sweep type)"))
save_image("HprimeFractionSigLargeEffect",PLOT)
