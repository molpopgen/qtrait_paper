library(dplyr)
library(readr)
library(tidyr)
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
    filter((g-5e4)/5e3<1e2/5e3,sweep_type != 'none',abs(s)>=2*sqrt(2)*sqrt(mu)) %>%
    group_by(mu,opt,repid,locus,sweep_type) %>% tally() %>%
    spread(sweep_type,n) %>% #This turns the sweep type into a column and assigns how many sweeps of each type to the locus
    right_join(dummy,by=c("mu","opt","locus","repid")) %>%
    mutate_all(funs(replace(.,is.na(.),0)))


sweed=read_delim("SweeD_Sig.txt.gz",delim=" ")
print(head(sweed))

power_by_sweep_type = sweed %>% left_join(fixations,by=c("mu","opt","repid","locus")) %>%
           mutate(sweep_type = case_when((hard > 0 & soft == 0) ~ 0,
                            (hard == 0 & soft > 0) ~ 1,
                            (hard > 0 & soft > 0) ~ 2,
                            (hard+soft == 0) ~ 3)) %>%
            group_by(mu,opt,generation,sweep_type) %>%
            summarise(psig_given_sweep = sum(sig05)/n()) %>%
            mutate(scaled_time = (generation-5e4)/5e3) %>%
            filter(scaled_time >= -0.4, scaled_time <= 1.0, sweep_type != 2)


save_image <- function(stat,img)
{
    trellis.device(device="pdf",file=paste(stat,".pdf",sep=""),height=10,width=10)
    trellis.par.set("fontsize",list(text=18))
    print(img)
    dev.off()
}
# KEY=list(space="top",columns=3,
#          title="Sweep category",
#          cex.title=1,#points=FALSE,
#          lines=list(lwd=rep(3,3)),
#          just=0.5,
#          text=list(c("New mutation","Standing var.","None")))
KEY=list(space="top",columns=3,title="Sweep category",
         cex.title=2,points=FALSE,just=0.5,lines=list(lwd=rep(3,3),lty=c(1,2,3)),
          text=list(c("New mutation","Standing var.","None")))
STRIP=strip.custom(strip.names = TRUE,sep=" = ", 
                   var.name = c(expression(mu),expression(z[o])),bg=c("white"))
bwtheme=standard.theme("pdf",color=FALSE)
bwtheme$fontsize$text = 18
Plot = xyplot(psig_given_sweep~scaled_time| as.factor(mu)*as.factor(opt), #:as.factor(sweep_type),
              group=sweep_type,lwd=2,
                  type='l',data=power_by_sweep_type,
                  par.settings=bwtheme,
                  key=KEY,
                  xlab="Time since optimum shift (units of N generations)",
                  ylab="Pr(significant|sweep category)",
                  scales=list(x=list(at=seq(-0.4,1.0,0.2),rot=45),alternating=F),
                  strip=STRIP,
                  xlim=c(-0.5,1)
                  )
save_image('FractionSweeDLargeEffect',Plot)
