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

    def getObjList(self, obj, exactQuery={}, approxQuery={}, orderBy="date"):
        """ Get a list of all objects satisfying the query """

        listRes = []

        # print(self.table.table.__dict__)

        # Convert bool to int
        queryStrings = []
        for key,value in exactQuery.items():
            if(value==None):
                queryStrings.append("({} IS NULL)".format(key))
            elif(type(value) == bool):
                queryStrings.append("({}={})".format(key,int(value)))
            else:
                queryStrings.append("({}='{}')".format(key,value))


        query = "SELECT {} FROM {} WHERE ".format(obj, self.table.table) +\
                " AND ".join(queryStrings)+" AND" +\
                " AND ".join(["({} BETWEEN {} AND {})".\
            format(key,value-tol,value+tol) for key,value in approxQuery.items()])+\
                " ORDER BY {} DESC".format(orderBy)

        for e in self.db.query(query):

            if obj=='eigv':
                listRes.append(scipy.fromstring(e[obj]).reshape(
                    e['neigs'], e['basisSize']))
            elif obj=='spec':
                listRes.append(sorted(scipy.fromstring(e[obj])))
            else:
                listRes.append(e[obj])

        return listRes

    def getEigs(self, k, ren, g, L, ET, EL=None, ELp=None, loc2=None, test=None,
            neigs=6, nonloc3mix=None, loc3mix=None, loc3=None):

        approxQuery = {"g":g, "L":L, "ET":ET}
        exactQuery = {"k": k, "ren":ren, "neigs":neigs, "finiteL":True}

        if ren=="rentails":
            if EL==None:
                EL = ratioELET*ET
            if ELp==None:
                ELp = ratioELpET*ET
            ELpp = ratioELppELp*ELp

            if loc2 != None:
                exactQuery["loc2"] = loc2
            if loc3mix != None:
                exactQuery["loc3mix"] = loc3mix
            if nonloc3mix != None:
                exactQuery["nonloc3mix"] = nonloc3mix
            if loc3 != None:
                exactQuery["loc3"] = loc3
            if test != None:
                exactQuery["test"] = test

            approxQuery["EL"] = EL
            approxQuery["ELp"] = ELp
            approxQuery["ELpp"] = ELpp

        try:
            ret = self.getObjList('spec', exactQuery, approxQuery)[0]

        except IndexError:
            msg = "L={}, ET={}, ren={}".format(approxQuery["L"],
                    approxQuery["ET"], exactQuery["ren"])
            print("Not found: ", exactQuery, approxQuery)
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
