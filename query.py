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
    values = ("k", "L", "ET", "g","neigs")

    # db = database.Database(dbname="spectraJson.db",useJson=True)
    db = database.Database()

    exactQuery = {}
    data = [db.getObjList(x, exactQuery=exactQuery, approxQuery=query) for x in values]
    res = set(zip(*data))

    for entry in sorted(res):
        print(", ".join(map(lambda x: ":".join((x[0],str(x[1]))), zip(values, entry))))

if __name__ == "__main__":
    main(sys.argv)
