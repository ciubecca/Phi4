import scipy
import json
import dataset
import datetime

dbname2="spectra2.db"
tablename2="spectra2"
db2 = dataset.connect('sqlite:///'+dbname2)
table2 = db2[tablename2]

dbname1="spectra.db"
tablename1="spectra"
db1 = dataset.connect('sqlite:///'+dbname1)
table1 = db1[tablename1]

for e in table1:
    e['eigv'] = scipy.array(json.loads(e['eigv'])).tostring()
    e['spec'] = scipy.array(json.loads(e['spec'])).tostring()
    table2.insert(e)
