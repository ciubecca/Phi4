import sys
import scipy
import math
from scipy import pi, log, array, sqrt
import database

dbname = "data/spectra3.db"

def main(argv):
    args = " <keyword:value> ..."

    query = {}
    queryStr = {"finiteL":True, "ren":"rentails"}

    for arg in argv[1:]:
        key, value = arg.split(':')
        try:
            query[key] = float(value)
        except ValueError:
            queryStr[key] = value

    # Hardcoded parameters
    values = ("k", "L", "ET", "g")

    # db = database.Database(dbname="spectraJson.db",useJson=True)
    db = database.Database(dbname)

    exactQuery = queryStr
    data = [db.getObjList(x, exactQuery=exactQuery, approxQuery=query) for x in values]
    res = set(zip(*data))

    for entry in sorted(res):
        print(", ".join(map(lambda x: ":".join((x[0],str(x[1]))), zip(values, entry))))

if __name__ == "__main__":
    main(sys.argv)
