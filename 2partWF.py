import sys
import matplotlib.pyplot as plt
import scipy
import math
import json
from scipy import pi, log, log10, array, sqrt, stats
from matplotlib import rc
import database
import finiteVolH

output = "svg"

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

def main(argv):
    args = " <g> <Emax> <L>"
    if len(argv) < 4:
        print(argv[0], args)
        return -1

    g = float(argv[1])
    Emax = float(argv[2])
    L = float(argv[3])

    params = {'legend.fontsize': 8}
    plt.rcParams.update(params)

    db = database.Database()

    exactQuery = {"ren":"raw", "k":1}
    approxQuery = {"g":g, "Emax":Emax, "L":L}
    eigv = db.getObjList("eigv", exactQuery=exactQuery, approxQuery=approxQuery)[0]

    vacuumEigv = eigv[0]
    print(vacuumEigv)

    # plt.figure(1)

    # data = vacuum["raw"]/xList
    # plt.plot(xList, data, linewidth=0.7,
            # dashes = [3,2], marker='.', markersize=3, color=evenColor, label=label("raw"))


#     plt.figure(1, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
    # #plt.xlim(min(xList)-0.01, max(xList)+0.01)
    # plt.title("m=1, dim0="+str(dim0)+", g="+str(g)+", n0eigs="+str(n0eigs)+", ordM="+str(ordMass))
    # plt.xlabel(r"$L$")
    # plt.ylabel(r"$\Lambda$")
    # plt.legend(loc='lower right', prop={'size':10})
    # plt.savefig("fig_vacVsL_g="+str(g)+"_n0eigs="+str(n0eigs)+"_dim0="+str(dim0)+"_ordM="+str(ordMass)+".svg")


if __name__ == "__main__":
    main(sys.argv)
