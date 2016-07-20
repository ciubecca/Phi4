import scipy
import json
import dataset
import datetime

dbnameJson = "spectraJson.db"
tablenameJson = "spectra"
dbJson = dataset.connect('sqlite:///'+dbnameJson)
tableJson = dbJson[tablenameJson]

dbname = "spectra.db"
tablename = "spectra"
db = dataset.connect('sqlite:///'+dbname)
table = db[tablename]

for e in tableJson:
    del e["id"]
    e['eigv'] = scipy.array(json.loads(e['eigv'])).tostring()
    e['spec'] = scipy.array(json.loads(e['spec'])).tostring()
    table.insert(e)
