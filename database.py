import scipy
import json
import dataset
import datetime

rentypes = ["raw","renlocal","rensubl"]

#FIXME Feature needed: forbid merging json and non-json data
class Database():
    def __init__(self, dbname="spectra.db", tablename="spectra", useJson=False):
        self.db = dataset.connect('sqlite:///'+dbname)
        self.table=self.db[tablename]
        self.useJson = useJson

    def insert(self, k, L, Emax, g, spec, eigv, basisSize, neigs, ren, cutoff=5.):
        if(basisSize*neigs != eigv.size):
            print(eigv.size)
            raise ValueError("basisSize, neigs and eigv dimension don't match")

        if ren not in rentypes:
            raise ValueError("ren argument must be in {}".format(", ".join(rentypes)))

        if self.useJson==True:
            self.table.insert(dict(date=datetime.datetime.now(), k=k, L=L, Emax=Emax, g=g, ren=ren, eigv=json.dumps(eigv.tolist()), \
                                cutoff=cutoff, spec=json.dumps(spec.tolist()), basisSize=basisSize, neigs=neigs))
        else:
            self.table.insert(dict(date=datetime.datetime.now(), k=k, L=L, Emax=Emax, g=g, ren=ren, eigv=eigv.tostring(), \
                                cutoff=cutoff, spec=spec.tostring(), basisSize=basisSize, neigs=neigs))

    # Get a list of all objects satisfying the query
    def getObjList(self, obj, exactQuery={}, approxQuery={}, boundQuery={}, orderBy=None):
        t1 = [e for e in self.table.find(**exactQuery)]
        listRes = []
        for e in t1:
            if all([abs(e[key]-value)<10.**(-12.) for key,value in approxQuery.items()]) and \
                all([value[0]<=e[key]<value[1] for key,value in boundQuery.items()]):
                if obj=='eigv':
                    if self.useJson == True:
                        listRes.append(json.loads(e[obj]))
                    else:
                        listRes.append(scipy.fromstring(e[obj]).reshape(e['neigs'], e['basisSize']))
                elif obj=='spec':
                    if self.useJson == True:
                        listRes.append(json.loads(e[obj]))
                    else:
                        listRes.append(scipy.fromstring(e[obj]))
                else:
                    listRes.append(e[obj])

        if orderBy==None:
            return listRes
        else:
            return [y for (x,y) in sorted(zip(orderBy, listRes))]


    def migrateJsonToBytes(self, otherdb)
        for e in self.table:
            del e["id"]
            e['eigv'] = scipy.array(json.loads(e['eigv'])).tostring()
            e['spec'] = scipy.array(json.loads(e['spec'])).tostring()
            otherdb.table.insert(e)
