library(dplyr)
library(tidyr)
library(DBI)
library(lattice)
library(viridis)

outer = c("main_results", "lowrec", "hirec")
rho = c(1e3,1e2,1e4)

data = tibble()
for (i in 1:length(outer))
{
        dbname = paste(outer[i],"_fixations.sqlite3",sep="")
        db = dbConnect(RSQLite::SQLite(), dbname)
        t = tbl(db, "data")

        df = collect(t) %>% 
            mutate(rho=rho[i])

        data = bind_rows(data, df)

        dbDisconnect(db)
}

data = data %>%
    mutate(rel_origin = g - 5e4,
           sojourn_time = ftime - g,
           scaled_gamma = 5000*s*s,
           strength = NA)

print(paste(min(data$s),max(data$s)))

data$strength[which(data$scaled_gamma >= 0 & data$scaled_gamma < 10)] = 1
data$strength[which(data$scaled_gamma >= 10 & data$scaled_gamma < 100)] = 2
data$strength[which(data$scaled_gamma >= 100 & data$scaled_gamma < 1000)] = 3
data$strength[which(data$scaled_gamma >= 1000)] = 4

data = data %>% filter(g + ftime >= 5e4, g <= 5e4+50)
# print(data)
# foo = data %>% filter(rho==1e4,mu==5e-3) 
# 
# print(foo[foo$s == min(foo$s),])
# print(foo[foo$repid==34,])

print(levels(factor(data$strength)))

COLORS=viridis(4)

STRIP=strip.custom(strip.names = TRUE,sep=" = ", 
                   var.name = c(expression(mu),expression(rho)), bg=c("white"))
KEYTEXT=list(c(expression(paste(0,"" <= "",'N',gamma^2,"" < "",10)),
expression(paste(10,"" <= "",'N', gamma^2,"" < "",100)),
 expression(paste(100,"" <= "",'N', gamma^2,"" < "",1000)), 
expression(paste('N',gamma^2,"" >= "",1000))))
KEY=list(space="top",columns=2,
         cex.title=0.25,
         points=list(pch=rep('.',length(COLORS)),col=COLORS,cex=rep(8,length(COLORS))),
         text=KEYTEXT)
         #text=list(sort(unique(data$strength))[c(4,1,2,3)]))
plt = xyplot(sojourn_time/5e3 ~ s|as.factor(mu)*as.factor(rho), data=data, 
             group=as.factor(strength),
             par.settings=simpleTheme(col=COLORS),
             pch='.',cex=5,
              scales=list(cex=1,alternating=F),
             xlab=expression(paste("Effect size(",gamma,")")),
             ylab="Fixation time (units of N generations)",
             ylim=c(0,2.6),
             key=KEY,
             strip=STRIP)

trellis.device(device='pdf',file='SojournTimesVaryRho.pdf',height=10,width=10)
trellis.par.set("fontsize",list(text=18))
print(plt)
dev.off()
