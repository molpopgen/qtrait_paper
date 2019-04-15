library(dplyr)
library(tidyr)
library(DBI)
library(ggplot2)
library(tibble)
library(viridis)

df = tibble()
rhovals = c(1e2,1e3,1e4)
I=1
twoN = 1e4
for(directory in c("lowrec", "main_results", "hirec"))
{
    db = dbConnect(RSQLite::SQLite(),paste(directory,"/himu/replicate1_ld.sqlite3",sep=""))
    dbt = tbl(db, 'data')
    for (minc2 in c(1,0.05*twoN, 0.25*twoN))
    {
        q = dbt %>% select(generation,pos1,pos2,c1,c2,D) %>%
            filter(c2 >= minc2, c2 <= twoN - minc2) %>%
            # mutate(intralocus = ifelse(abs(pos1-pos2) < 1.0, 1, 0)) %>%
            # group_by(generation, pos1, intralocus) %>% summarise(mD=mean(D))
            group_by(generation, pos1) %>% summarise(mD=mean(D))
        t = collect(q) %>% mutate(rho = rhovals[I],minc2 = minc2/twoN)
        df = bind_rows(df, t)
    }
    dbDisconnect(db)
    I <- I + 1
}

label_rho <- function(x)
{
    return (expression(paste(rho," = ",x)))
}

label_count <- function(x)
{
    return (paste("MAF =",as.numeric(x)/twoN,sep=' '))
}

p = ggplot(df,aes(x=generation-5e4,y=mD,color=as.factor(pos1),alpha=0.5)) + geom_line(alpha=0.5) + 
    #facet_grid(minc2~as.factor(rho),labeller=labeller(minc2=label_count,rho=as_labeller(label_rho))) + 
    facet_grid(minc2~rho,
               labeller=label_bquote(cols=rho ~"="~ .(rho),rows="MAF = " ~ .(minc2))) +
    theme_bw(base_size=14) +
    theme(legend.position='none',axis.text.x=element_text(angle=45,vjust=1,hjust=1)) + 
    xlab("Generations since optimum shift") + 
    ylab("Mean disequilibrium statistic, D")

ggsave("LDsinglerep.pdf",height=10,width=10,units="in")

# up = unique(df$pos1)
# 
# t = df %>% filter(pos1==up[1])
# plot(t$generation, t$mD)
# for (i in up[2:length(up)])
# {
#     t = df %>% filter(pos1==i)
#     lines(t$generation, t$mD)
# }
# 
