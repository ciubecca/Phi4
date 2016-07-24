import sys
import scipy
import math
import json
from scipy import pi, log, array, sqrt
from matplotlib import rc
import database
import finiteVolH


def main(argv):
    args = " <keyword:value> ..."

    query = {}

    for arg in argv[1:]:
        key, value = arg.split(':')
        query[key] = float(value)

    # Hardcoded parameters
    values = ("L", "Emax", "basisSize", "g", "occmax")

    # db = database.Database(dbname="spectraJson.db",useJson=True)
    db = database.Database()

    exactQuery = {}
    data = [db.getObjList(x, exactQuery=exactQuery, approxQuery=query) for x in values]
    # Check integrity of eigv data
    eigv = db.getObjList("eigv", exactQuery=exactQuery, approxQuery=query)
    res = set(zip(*data))
    # Transpose list and remove duplicates
    # res = zip(*data)

    for x in sorted(res):
        print(list(zip(values, x)))

if __name__ == "__main__":
    main(sys.argv)
