import sys
import scipy
import math
from scipy import pi, log, array, sqrt
import database

ETmax = {5:32, 5.5:30, 6:28, 6.5:26.5, 7:25, 7.5:24, 8:23,
                8.5:22, 9:21, 9.5:20.5, 10:20}
ETmin = {5:10, 5.5:10, 6:10, 6.5:10, 7:10, 7.5:10, 8:10,
                8.5:10, 9:10, 9.5:10, 10:10}

Llist = sorted(ETmax.keys())

dbname = "data/spectra3.db"
dbname = "data/spectra.db"

def main(argv):

    query = {}
    exactQuery = {"finiteL":True, "ren":"rentails"}

    db = database.Database(dbname)

    if len(argv)<2:
        print(argv[0], " <g>")
        sys.exit(1)

    g = float(argv[1])

    for L in Llist:
        for ET in scipy.linspace(ETmin[L],ETmax[L],(ETmax[L]-ETmin[L])*2+1):
            for k in (-1,1):
                exactQuery["k"] = k
                approxQuery = {"L":L, "ET":ET, "g":g}
                try:
                    db.getObjList("spec",exactQuery=exactQuery,
                            approxQuery=approxQuery)[0]
                except IndexError:
                    print("k={}, L={}, ET:{}, g={} not found".format(k,L,ET,g))


if __name__ == "__main__":
    main(sys.argv)
