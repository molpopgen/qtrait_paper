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
    outfile = infile.replace(".db","_recent_fixations.db")

    if os.path.exists(outfile):
        os.remove(outfile)

    oconn = sqlite3.connect(outfile)
    x = infile.split("_")
    opt = float(x[3])
    mu = float(x[5])
    conn = sqlite3.connect(infile)
    c = conn.cursor()
    for i in range(256):
        c.execute(
            'create temp table mfreq as select *,MAX(freq) as mfreq,MAX(generation) as mgen from freqs where repid == {0} group by repid,origin,pos,esize'.format(i))
        c.execute("create temp table fix as select * from mfreq where mfreq == 1")
        c.execute("create temp table fixtraj as select * from freqs inner join fix on fix.pos == freqs.pos and fix.esize == freqs.esize and fix.origin == freqs.origin and fix.repid == freqs.repid and fix.mgen >= 40000 and fix.origin <= 75000")
        conn.commit()

        x = pd.read_sql(
            'select generation,pos,esize,origin,freq,repid from fixtraj where repid == {0}'.format(i), conn)
        x['opt'] = [opt] * len(x.index)
        x['mu'] = [mu] * len(x.index)
        x.to_sql('freqs', oconn, if_exists='append')
        c.execute('drop table mfreq')
        c.execute('drop table fix')
        c.execute('drop table fixtraj')
        conn.commit()
    conn.close()
    oconn.close()
    return None

if __name__ == "__main__":
    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = {executor.submit(get_recent_fixations, i): i for i in sys.argv[1:]}
        for fut in concurrent.futures.as_completed(futures):
            fn = fut.result()
