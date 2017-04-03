import database
import extrapolate

db1 = database.Database("data/spectra.db")
db2 = database.Database("data/spectra2.db")
db3 = database.Database("data/spectra3.db")

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
