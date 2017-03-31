import scipy
import json
import dataset
import datetime
from scipy import array
import sys

rentypes = ["raw","renloc","rentails"]

# Ratio between EL and ET
ratioELET = 3
# Ratio between ELp and ET
ratioELpET = 2
# Ratio between ELpp and ELp
ratioELppELp = 1.5

# The range where to look for numerical values in the database
tol = 10**-6

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
    def getObjList(self, obj, exactQuery={}, approxQuery={}, orderBy="date"):

        listRes = []

        # Convert bool to int
        exactQueryNoBool = {}
        for key,value in exactQuery.items():
            if(type(value) == bool):
                exactQueryNoBool[key] = int(value)
            else:
                exactQueryNoBool[key] = value

        exactQueryStr = " AND ".join("({}='{}')".\
                format(key,value) for key,value in exactQueryNoBool.items())

        # print(self.table.table)
        query = "SELECT {} FROM {} WHERE {} AND ".\
                format(obj, self.table.table, exactQueryStr)+\
                " AND ".join("({} BETWEEN {} AND {})".\
                format(key,value-tol,value+tol) for key,value in approxQuery.items())+\
                " ORDER BY {}".format(orderBy)
        # print(query)
        result = self.db.query(query)

        # for e in self.table.find(**exactQuery):
        for e in result:

            if obj=='eigv':
                listRes.append(scipy.fromstring(e[obj]).reshape(
                    e['neigs'], e['basisSize']))
            elif obj=='spec':
                listRes.append(sorted(scipy.fromstring(e[obj])))
            else:
                listRes.append(e[obj])

        return listRes

    def getEigs(self, k, ren, g, L, ET, neigs=6):


        approxQuery = {"g":g, "L":L, "ET":ET}
        exactQuery = {"k": k, "ren":ren, "neigs":neigs, "finiteL":True}

        if ren=="rentails":
            EL = ratioELET*ET
            ELp = ratioELpET*ET
            ELpp = ratioELppELp*ELp


            exactQuery["maxntails"] = None
            # exactQuery["tailsComputedAtET"] = ET
            approxQuery["EL"] = EL
            approxQuery["ELp"] = ELp
            approxQuery["ELpp"] = ELpp

        try:
            ret = self.getObjList('spec', exactQuery, approxQuery)[0]

        except IndexError:
            print("Not found:", exactQuery, approxQuery)
            exit(-1)
        except TypeError as e:
            print(exactQuery, approxQuery)
            raise e

        return ret

    # Convert to json format
    def convert(self, newdbname):
        newdb = Database(newdbname)
        for e in self.table:
            e["spec"] = json.dumps(sorted(scipy.fromstring(e["spec"])))
            newdb.table.insert(e)
