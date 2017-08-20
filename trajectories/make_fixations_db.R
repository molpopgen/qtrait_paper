library(data.table)
library(dplyr)
library(RSQLite)

n=commandArgs(trailing=TRUE)
fixations_db = src_sqlite("output/fixations.db",create=TRUE)

for(ni in n)
{
    x=fread(paste("gunzip -c",ni,sep=" "))
    dbWriteTable(fixations_db$con,"fixations",as.data.frame(x),append=TRUE)
}
