# Defines class used to store computed eigenvalues and eigenvectors

import scipy
import json
import dataset
import datetime
from scipy import array
import sys

# Tolerance parameter when looking for numerical values in the database
tol = 10**-6

class Database():
    def __init__(self, dbname="data/spectra.db", tablename="spectra"):
        self.db = dataset.connect('sqlite:///'+dbname)
        self.table = self.db[tablename]

    def insert(self, datadict):
        """ Insert entry in the database
        datadict: dictionary object, in the form {column: value}
        """
        datadict["date"] = datetime.datetime.now()
        datadict["spec"] = datadict["spec"].tostring()
        self.table.insert(datadict)

    def getObjList(self, obj, exactQuery={}, approxQuery={}, orderBy="date"):
        """ Get a list of all objects satisfying the query
        obj: which column of the database to extract
        exactQuery: query in the form {column: value}, where value must match exactly
        approxQuery: query in the form {column: value}, where value must match approximately
        """

        listRes = []

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
                listRes.append(scipy.fromstring(e[obj]))

            elif obj=='spec':
                listRes.append(sorted(scipy.fromstring(e[obj])))
            else:
                listRes.append(e[obj])

        return listRes


    def getEigs(self, k, ren, g2, g4, L, ET, Lambda, neigs=6):
        """ Extract eigenvalues from the database corresponding to given parameters
        k: parity quantum number
        ren: type of renormalization
        g2: phi^2 coupling
        g4: phi^4 coupling
        L: torus side
        ET: energy cutoff
        Lambda: momentum cutoff
        neigs: number of eigenvalues computed during diagonalization
        """

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

    def convert(self, newdbname):
        """ Convert database to json format """
        newdb = Database(newdbname)
        for e in self.table:
            e["spec"] = json.dumps(sorted(scipy.fromstring(e["spec"])))
            newdb.table.insert(e)
