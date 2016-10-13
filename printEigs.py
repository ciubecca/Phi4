import sys
import scipy
import math
from scipy import pi, log, log10, array, sqrt, stats
import database

renlist = ("raw", "renloc", "rentails")

neigs = 10

def EL(ET):
    return ET*2.5

def plotvsE(Elist):

    db = database.Database()

    exactQuery = {}
    approxQuery = {"g":g, "L":L}

    E0 = {}
    E1 = {}
    for ren in renlist:
        E0[ren] = []
        E1[ren] = []

        for ET in ETlist:
            exactQuery["ren"] = ren
            approxQuery["ET"] = ET
            exactQuery["k"] = 1
            if ren=="renloc":
                approxQuery["EL"] = ET
            elif ren=="rentails":
                approxQuery["EL"] = EL(ET)
                # approxQuery["EL"] = max(ETlist)+ELETdiff
            E0[ren].append(db.getObjList('spec', exactQuery, approxQuery)[0][0])

            exactQuery["k"] = -1
            E1[ren].append(db.getObjList('spec', exactQuery, approxQuery)[0][0])


    # VACUUM ENERGY
    data = E0["raw"]
    print("raw vacuum energy")
    print(data)

    data = E0["renloc"]
    print("renloc vacuum energy")
    print(data)

    data = E0["rentails"]
    print("rentails vacuum energy")
    print(data)

    # MASS
    data = array(E1["raw"])-array(E0["raw"])
    print("raw mass")
    print(data)

    data = array(E1["renloc"])-array(E0["renloc"])
    print("renloc mass")
    print(data)

    data = array(E1["rentails"])-array(E0["rentails"])
    print("rentails mass")
    print(data)



argv = sys.argv
if len(argv) < 5:
    print(argv[0], "<L> <g> <ETmin> <ETmax>")
    sys.exit(-1)

L = float(argv[1])
g = float(argv[2])
ETmin = float(argv[3])
ETmax = float(argv[4])

ETlist = scipy.linspace(ETmin, ETmax, 2*(ETmax-ETmin)+1)
print("ETlist:", ETlist)



plotvsE(ETlist)

