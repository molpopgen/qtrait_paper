import sqlite3
import pandas
import glob

for f in glob.glob('*.sqlite3'):
    with sqlite3.connect(f) as conn:
        df = pandas.read_sql('select max(abs(esize)) as esize from data', conn)
        print(f, df.esize.iloc[0])

