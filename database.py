import scipy
import json
import dataset
import datetime

class Database():
    def __init__(self, dbname="spectra2.db", tablename="spectra2"):
        self.db = dataset.connect('sqlite:///'+dbname)
        self.table=self.db[tablename]

    def insert(self, k, L, Emax, g, spec, eigv, basisSize, neigs, ren, cutoff=5.):

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
                    listRes.append(scipy.fromstring(e[obj]).reshape(e['neigs'],-1))
                elif obj=='spec':
                    listRes.append(scipy.fromstring(e[obj]))
                else:
                    listRes.append(e[obj])

        if orderBy==None:
            return listRes
        else:
            return [y for (x,y) in sorted(zip(orderBy, listRes))]


