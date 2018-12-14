library(dplyr)
library(stringr)
library(viridis)
library(tidyr)
library(lattice)

process_files <- function(d, param)
{
    df=NA
    for (i in d)
    {
        print(i)
        x = str_split(i,regex('\\.'))[[1]]
        mu = as.numeric(str_replace(paste(x[2],x[3],sep="."),'mu',''))
        plarge = as.numeric(str_replace(paste(x[4],x[5],sep="."),'plarge',''))
        ghat = 2*sqrt(2)*sqrt(mu)
        db = DBI::dbConnect(RSQLite::SQLite(),i)
        dbt=tbl(db,'data')
        results = collect(dbt) %>%
            group_by(generation) %>%
            summarise(mz = mean(zbar),
                      mvg = mean(vg))
        results=results%>%
            mutate(mu=mu,plarge=plarge,des=param)
        if(is.na(df))
        {
            df=results
        } else
        {
            df=rbind(df,results)
        }
        DBI::dbDisconnect(db)
    }
    df=df%>% mutate(scaled_time=(generation-5e4)/5e3)
    df
}
gamma_small_files = dir('.',pattern='qtrait_gamma_small\\..+sqlite3')
gamma_files = dir('.',pattern='qtrait_gamma\\..+sqlite3')
gauss_files = dir('.',pattern='qtrait\\..+sqlite3')


# NOTE: we must be careful about labels across all scripts...
gauss_df = process_files(gauss_files, 'gauss')
gamma_df = process_files(gamma_files, 'gamma1')
gamma2_df = process_files(gamma_small_files, 'gamma2')

df = rbind(gauss_df, gamma_df, gamma2_df) 

COLORS=rev(viridis(3))
STRIP=strip.custom(strip.names = TRUE,sep=" = ", 
                   var.name = c(expression(mu),expression(paste("Pr(|",gamma,"| >= ",hat(gamma),")"))),bg=c("white"))
key_text=list(c(expression(paste(Gamma,", shape = 1.0")),
              expression(paste(Gamma,", shape = 0.5")),
              expression(paste("gaussian"))))
KEY=list(space="top",columns=3,
         title="Distribution of effect sizes",
         cex.title=1,#points=FALSE,
         lines=list(lwd=rep(3,length(COLORS)),col=COLORS),
         just=0.5,
         text=key_text)
p = xyplot(mz~scaled_time|as.factor(mu)*as.factor(plarge),
           data=df,
           type='l',
           par.settings=simpleTheme(col=COLORS),
           scales=list(cex=1,alternating=F,x=list(rot=45)),
           strip=STRIP,
           group=des,
           lwd=3,
           xlim=c(-0.1,0.2),
           ylim=c(0.,1.1),
           xlab="Time since optimum shift (units of N generations)",
           ylab="Mean genetic value",
           key=KEY)
trellis.device('pdf',file='MeanGeneticValueVaryPlarge.pdf',height=10,width=10)
trellis.par.set("fontsize",list(text=18))
print(p)
dev.off()

p = xyplot(log10(mvg)~scaled_time|as.factor(mu)*as.factor(plarge),
           data=df,
           type='l',
           par.settings=simpleTheme(col=COLORS),
           scales=list(cex=1,alternating=F,x=list(rot=45)),
           strip=STRIP,
           group=des,
           lwd=3,
           xlim=c(-0.1,0.2),
           #ylim=c(0.,0.6),
           xlab="Time since optimum shift (units of N generations)",
           ylab="log10(Mean genetic variance)",
           key=KEY)
trellis.device('pdf',file='MeanGeneticVarianceVaryPlarge.pdf',height=10,width=10)
trellis.par.set("fontsize",list(text=18))
print(p)
dev.off()

