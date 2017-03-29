import database

db1 = database.Database("data/spectra.db")
db2 = database.Database("data/spectra2.db")
db3 = database.Database("data/spectra3.db")

for x in db1.table:
    del x['id']
    db3.table.insert(x)

for x in db2.table:
    del x['id']
    db3.table.insert(x)
