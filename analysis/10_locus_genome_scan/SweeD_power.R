library(dplyr)
library(tidyr)
library(readr)
library(viridis)
library(lattice)

x = read_delim("raw_sweed.txt.gz",delim=" ")

db = src_sqlite("sweed_nulldist.db")
dt = tbl(db,'data')

ndist = collect(dt)

cvals = as.numeric(quantile(ndist$LR,prob=c(0.95,0.99,1-1e-4)))
x = x %>% mutate(sig05 = ifelse(LR >= cvals[1],1,0),
                 sig01 = ifelse(LR >= cvals[2],1,0),
                 sig00001 = ifelse(LR >= cvals[3],1,0))

# Write the intermediate file, as this could be useful later
f = gzfile("SweeD_Sig.txt.gz","wb")
write_delim(x,f,delim=" ")
close(f)

# Define power as ">= 1 significant locus per replicate by generation"

print(head(x))
# power = x %>%
#  group_by(generation,repid,mu,opt) %>%
#  summarise(sum_sig05 = ifelse(sum(sig05) > 0, 1, 0),
#         sum_sig01 = ifelse(sum(sig01) > 0, 1, 0),
#         sum_sig00001 = ifelse(sum(sig00001) > 0, 1, 0))
# 
power = x %>%
    group_by(generation,mu,opt) %>% 
    summarise(power05 = sum(sig05)/n(),
           power01 = sum(sig01)/n(),
           power00001 = sum(sig00001)/n()) %>%
    gather(alpha_rej,power,c(power05,power01,power00001)) %>%
    mutate(alpha_rej=as.numeric(sub('power0', '0.0', alpha_rej)),
           scaled_time=(generation-50000)/50000)

COLORS=viridis(length(cvals))
KEY=list(space="top",columns=3,title="Nominal significance level",
         cex.title=1,points=FALSE,lines=TRUE,just=0.5)
STRIP=strip.custom(strip.names = TRUE, 
                   var.name = c(expression(mu),expression(z[o])),bg=c("white"))
powerPlot = xyplot(power ~ scaled_time| as.factor(mu)*as.factor(opt),group=alpha_rej,
                  type='l',data=power,
                  par.settings=simpleTheme(col=COLORS),
                  auto.key=KEY,
                  xlab="Time since optimum shift (units of N generations)",
                  ylab="Power of composite likelihood ratio test",
                  scales=list(cex=1,alternating=F),
                  strip=STRIP)

trellis.device(device="pdf",file=paste("SweeDRawPower.pdf",sep=""),height=10,width=10)
print(powerPlot)
dev.off()
