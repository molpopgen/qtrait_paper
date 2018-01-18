import sqlite3
import sys
import os
import pandas as pd
import concurrent.futures

# Process the big trajectory databasees
# and pull out just the subset for fixations
# occurring more recently into the simulation.

def get_recent_fixations(args):
    infile = args


    x = infile.split("_")
    opt = float(x[3])
    mu = float(x[5])
    conn = sqlite3.connect(infile)
    c = conn.cursor()
    rv = []
    for i in range(256):
        c.execute(
            'create temp table mfreq as select *,MAX(freq) as mfreq,MAX(generation) as mgen from freqs where repid == {0} group by repid,origin,pos,esize'.format(i))
        c.execute("create temp table fix as select * from mfreq where mfreq == 1")
        c.execute("create temp table fixtraj as select * from freqs inner join fix on fix.pos == freqs.pos and fix.esize == freqs.esize and fix.origin == freqs.origin and fix.repid == freqs.repid ")
        conn.commit()

        x = pd.read_sql(
            'select generation,pos,esize,origin,freq,repid from fixtraj where repid == {0}'.format(i), conn)
        x = x[['origin','pos','esize','repid']].groupby(['origin','pos','esize']).agg(['count'])
        x.columns = x.columns.droplevel(0)
        x.rename(columns={'count':'sojourn_time'},inplace=True)
        x['opt'] = [opt] * len(x.index)
        x['mu'] = [mu] * len(x.index)
        x['repid'] = [i]*len(x.index)
        x.reset_index(inplace=True)
        rv.append(x)
        c.execute('drop table mfreq')
        c.execute('drop table fix')
        c.execute('drop table fixtraj')
        conn.commit()
    conn.close()
    return rv

if __name__ == "__main__":
    outfile = "all_fixation_summaries.sqlite3"
    if os.path.exists(outfile):
        os.remove(outfile)
    conn = sqlite3.connect(outfile)
    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = {executor.submit(get_recent_fixations, i): i for i in sys.argv[1:]}
        for fut in concurrent.futures.as_completed(futures):
            fn = fut.result()
            df = pd.concat(fn)
            df.to_sql('fixations',conn,index=False,if_exists='append')
    conn.close()
