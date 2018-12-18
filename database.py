import scipy
import json
import dataset
import datetime
from scipy import array
import sys

# The range where to look for numerical values in the database
tol = 10**-6

class Database():
    def __init__(self, dbname="data/spectra.db", tablename="spectra"):
        self.db = dataset.connect('sqlite:///'+dbname)
        self.table = self.db[tablename]

    def insert(self, datadict):

        datadict["date"] = datetime.datetime.now()
        datadict["spec"] = datadict["spec"].tostring()
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

            if obj=='eigv' and e[obj]!=None:
                # listRes.append(scipy.fromstring(e[obj]).reshape(
                    # e['neigs'], e['basisSize']))
# XXX Only vacuum state
                listRes.append(scipy.fromstring(e[obj]))

            elif obj=='spec':
                listRes.append(sorted(scipy.fromstring(e[obj])))
            else:
                listRes.append(e[obj])

        return listRes

    def getEigs(self, k, ren, g2, g4, L, ET, Lambda, neigs=6):

        approxQuery = {"g4":g4, "g2":g2, "L":L, "ET":ET, "Lambda":Lambda}
        exactQuery = {"k": k, "ren":ren, "neigs":neigs}

        try:
            ret = self.getObjList('spec', exactQuery, approxQuery)[0]

        except IndexError:
            msg = "L={}, ET={}, Lambda={}".format(approxQuery["L"],
                    approxQuery["ET"], exactQuery["Lambda"])
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
