import sys
import scipy
import math
from scipy import pi, log, array, sqrt
import database


def main(argv):
    args = " <keyword:value> ..."

    query = {}

    for arg in argv[1:]:
        key, value = arg.split(':')
        query[key] = float(value)

    # Hardcoded parameters
    values = ("L", "ET", "basisSize", "g", "EL", "spec")

    # db = database.Database(dbname="spectraJson.db",useJson=True)
    db = database.Database()

    exactQuery = {}
    data = [db.getObjList(x, exactQuery=exactQuery, approxQuery=query) for x in values]
    # Check integrity of eigv data
    # spec = db.getObjList("spec", exactQuery=exactQuery, approxQuery=query)
    # print(spec[0])
    res = set(zip(*data))
    # Transpose list and remove duplicates
    # res = zip(*data)

    for x in sorted(res):
        print(list(zip(values, x)))

if __name__ == "__main__":
    main(sys.argv)
