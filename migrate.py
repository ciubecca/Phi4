import database

db = database.Database(dbname="spectraJson.db", useJson=True)
otherdb = database.Database(dbname="spectra.db")

db.migrateJsonToBytes(otherdb)
