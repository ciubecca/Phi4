import scipy
import json
import dataset
import datetime
from scipy import array

rentypes = ["raw","renloc","rentails"]

# Ratio between EL and ET
ratioELET = 3
# Ratio between ELp and ET
ratioELpET = 2
# Ratio between ELpp and ELp
ratioELppELp = 1.5



#FIXME Feature needed: forbid merging json and non-json data
class Database():
    def __init__(self, dbname="data/spectra.db", tablename="spectra"):
        self.db = dataset.connect('sqlite:///'+dbname)
        self.table = self.db[tablename]

    def insert(self, datadict, spec):

        datadict["date"] = datetime.datetime.now()
        datadict["spec"] = spec.tostring()
        self.table.insert(datadict)

    # Get a list of all objects satisfying the query
    def getObjList(self, obj, exactQuery={}, approxQuery={}, boundQuery={},
            orderBy=None):

        listRes = []
        orderKey = []

        t1 = [e for e in self.table.find(**exactQuery)]

        for e in t1:
            if all([abs(e[key]-value)<10.**(-12.)
                    for key,value in approxQuery.items()]) \
                and all([value[0]<=e[key]<=value[1]
                    for key,value in boundQuery.items()]):


                if obj=='eigv':
                    listRes.append(scipy.fromstring(e[obj]).reshape(
                        e['neigs'], e['basisSize']))
                elif obj=='spec':
                    listRes.append(sorted(scipy.fromstring(e[obj])))
                else:
                    listRes.append(e[obj])

                if orderBy!=None:
                    orderKey.append(e[orderBy])

        if orderBy==None:
            return listRes
        else:
            return [y for (x,y) in sorted(zip(orderKey, listRes),
                key=lambda pair:pair[0])]


    def getEigs(self, k, ren, g, L, ET, neigs=6):

        EL = ratioELET*ET
        ELp = ratioELpET*ET
        ELpp = ratioELppELp*ELp

        approxQuery = {"g":g, "L":L, "ET":ET}
        exactQuery = {"k": k, "ren":ren, "neigs":neigs}
        boundQuery = {}

        if ren=="rentails":
            exactQuery["maxntails"] = None
            exactQuery["tailsComputedAtET"] = ET
            approxQuery["EL"] = EL
            approxQuery["ELp"] = ELp
            approxQuery["ELpp"] = ELpp

        try:
            ret = self.getObjList('spec', exactQuery, approxQuery, boundQuery,
                    orderBy="date")[0]

        except IndexError:
            print("Not found:", exactQuery, approxQuery)
            exit(-1)

        return ret
