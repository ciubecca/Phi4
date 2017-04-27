import sqlite3
import sys
import extrapolate

db1 = sqlite3.connect("data/spectra.db")
db2 = sqlite3.connect("data/spectra2.db")
db3 = sqlite3.connect("data/spectra3.db")

cursor1 = db1.cursor()
cursor1.execute('SELECT * FROM spectra')
output1 = cursor1.fetchall()

cursor2 = db2.cursor()
cursor2.execute("SELECT * FROM spectra WHERE (ren='renlocal' OR ren='raw')")
output2 = cursor2.fetchall()

print([description[0] for description in cursor1.description])
print(output1[0])
print([description[0] for description in cursor2.description])
print(output2[0])


sys.exit(0)

for x in db1.table:
    del x['id']
    db3.table.insert(x)

for x in db2.table:
# Don't trust rentails data
    if x['ren'] == 'rentails':
        continue
# Trust only data with high Emax
    if x['ET'] <= extrapolate.ETmax['rentails'][x['L']]:
        continue
    del x['id']
    db3.table.insert(x)
