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
        q = dbt %>% select(generation,pos1,pos2,e1,e2,c1,c2,D) %>%
            filter(generation <= 60000, c2 >= minc2, c2 <= twoN - minc2, 0.5*twoN*e1*e1 >= 10) %>%
            mutate(intralocus = ifelse(abs(pos1-pos2) < 1.0, 1, 0),
                   positive = ifelse(e1 > 0, '+', '-'),
                   positive2 = ifelse(e2 > 0, '+', '-')) %>%
            group_by(generation, pos1, positive, positive2, intralocus) %>% summarise(mD=mean(D),mine=min(e2),maxe=max(e2))
            #group_by(generation, pos1, positive) %>% summarise(mD=mean(D))
            #group_by(generation, pos1) %>% summarise(mD=mean(D))
        t = collect(q) %>% mutate(rho = rhovals[I],minc2 = minc2/twoN)
        #print(head(t))
        # x = t %>% group_by(pos1, positive, positive2, intralocus) %>% summarise(n=n())
        # print(x)
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

df = df %>% filter(intralocus == 1)
p = ggplot(df,aes(x=generation-5e4,y=mD,color=as.factor(pos1),shape=positive2)) + geom_point(alpha=0.25,size=0.05) + 
    scale_color_viridis(discrete=TRUE) +
    #facet_grid(minc2~as.factor(rho),labeller=labeller(minc2=label_count,rho=as_labeller(label_rho))) + 
    facet_grid(minc2*positive2~rho*positive,
               labeller=label_bquote(cols=rho ~"="~ .(rho)~","~.(positive),rows="MAF >="~.(minc2) ~", "~.(positive2))) +
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
