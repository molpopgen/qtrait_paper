import sqlite3
import sys

dbfile = sys.argv[1]

conn = sqlite3.connect(dbfile)

c = conn.cursor()

#Some DB setting for speed
c.execute("PRAGMA count_changes=OFF")
c.execute("PRAGMA journal_mode=MEMORY")
c.execute("PRAGMA cache_size(10000)")
#c.execute("PRAGMA temp_store=MEMORY")

#Create the table indexes
#CREATE TABLE freqs(generation int NOT NULL, origin int NOT NULL, pos real NOT NULL, esize real NOT NULL, freq real NOT NULL, repid int)
c.execute("create index if not exists repid on freqs (repid)")
c.execute("create index if not exists repid_origin on freqs (repid,origin)")
