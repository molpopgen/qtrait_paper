library(dplyr)

patterns = c("*.genome_scan.db","*.vg_over_time","*_sweed.db","*.genetic_values_per_locus.db")

disappointments = array()
I=1
for (p in patterns)
{
    files = dir('.',pattern=p)
    for (f in files)
    {
        #db = src_sqlite(f)
        db = DBI::dbConnect(RSQLite::SQLite(),f)
        dbt = tbl(db,'data')

        mgen = dbt %>%
            group_by(repid) %>%
            summarise(mg = max(generation))
        df = collect(mgen)
        mg = sort(unique(df$mg))
        if(mg[1] < 1e5)
        {
            disappointments[I]=f
            I=I+1
        }
        DBI::dbDisconnect(db)
    }
}
print(sort(disappointments))
