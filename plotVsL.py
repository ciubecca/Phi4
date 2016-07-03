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
version = "v3_5-LV"
eigTypes = ("raw", "loc")

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

def mkink(g):
    return 1/(12.*g) - (3./(2.*pi)-1./(4.*sqrt(3.)))

def main(argv):
    args = " <g> <n0eigs> <ordMass> <dim0>"
    if len(argv) < 5:
        print argv[0], args
        return -1

    g = float(argv[1])
    n0eigs = int(argv[2])
    ordMass = float(argv[3])
    dim0 = int(argv[4])

    # Hardcoded parameters
    neigs = 4
    xList =scipy.linspace(5,25,21)
    print xList

    plt.figure(1) # E_0
    plt.clf()

    plt.figure(2) # E_i-E_0
    plt.clf()

    params = {'legend.fontsize': 8}
    plt.rcParams.update(params)

    db = database.Database()

    oddColor = 'green'
    evenColor = 'purple'

    exactQuery = {"ren":"raw", "k":1, "n0eigs":n0eigs, "dim0":dim0, "finiteVolCouplings":True, "version":version, "ordMass":ordMass}
    approxQuery = {"g":g}
    boundQuery = {"neigs":[neigs,1000], "basisSize":[5000,10**5]}

    eigenvalues = {}
    for ren in eigTypes:
        eigenvalues[ren] = {}
        for k in (1,-1):
            eigenvalues[ren][k] = []
            for L in xList:
                exactQuery = {"ren":ren, "k":k, "n0eigs":n0eigs, "dim0":dim0, "finiteVolCouplings":True, "version":version, "ordMass":ordMass}
                approxQuery["L"] = L
                eigenvalues[ren][k].append(db.getObj('spec', exactQuery, approxQuery, boundQuery, maxAttr="basisSize"))


    vacuum = { }
    evenSpectrum = { }
    oddSpectrum = { }
    for ren in eigTypes:
        vacuum[ren] = scipy.array([e[0] for e in eigenvalues[ren][1]])
        evenSpectrum[ren] = [ scipy.array([e[i] for e in eigenvalues[ren][1]]) - vacuum[ren] for i in range(1,neigs) ]
        oddSpectrum[ren] = [ scipy.array([e[i] for e in eigenvalues[ren][-1]]) - vacuum[ren] for i in range(neigs) ]


    def label(ren, k=None):
        p = ""
        if k==1:
            p = r"$, Z_2=+$"
        elif k==-1:
            p = r"$, Z_2=-$"
        return ren+p


    # VACUUM ENERGY DENSITY
    plt.figure(1)

    # RAW
    data = vacuum["raw"]/xList
    plt.plot(xList, data, linewidth=0.7,
            dashes = [3,2], marker='.', markersize=3, color=evenColor, label=label("raw"))

    # LOC
    data = vacuum["loc"]/xList
    plt.plot(xList, data, linewidth=1.,
            dashes = [4,1], marker='+', markersize=3, color=evenColor, label=label("loc"))


    # MASS

    # RAW
    plt.figure(2)
    data = oddSpectrum["raw"][1]
    l=label("raw", -1)
    plt.plot(xList, data, linewidth=0.7,
            dashes = [3,2], marker='.', markersize=3, color=oddColor, label=l)
    data = evenSpectrum["raw"][0]
    l=label("raw", 1)
    plt.plot(xList, data, linewidth=0.7,
            dashes = [3,2], marker='.', markersize=3, color=evenColor, label=l)

    # LOC
    data = oddSpectrum["loc"][1]
    l=label("loc", -1)
    plt.plot(xList, data, linewidth=0.7,
            dashes = [4,1], marker='+', markersize=3, color=oddColor, label=l)

    data = evenSpectrum["loc"][0]
    l=label("loc", 1)
    plt.plot(xList, data, linewidth=0.7,
            dashes = [4,1], marker='+', markersize=3, color=evenColor, label=l)


    # E_2 - 2 MASS

    # RAW
    plt.figure(3)
    data = oddSpectrum["raw"][2]-2*oddSpectrum["raw"][1]
    l=label("raw", -1)
    plt.plot(xList, data, linewidth=0.7,
            dashes = [3,2], marker='.', markersize=3, color=oddColor, label=l)

    data = evenSpectrum["raw"][1]-2*evenSpectrum["raw"][0]
    l=label("raw", 1)
    plt.plot(xList, data, linewidth=0.7,
            dashes = [3,2], marker='.', markersize=3, color=evenColor, label=l)

    # LOC
    data = oddSpectrum["loc"][2]-2*oddSpectrum["loc"][1]
    l=label("loc", -1)
    plt.plot(xList, data, linewidth=0.7,
            dashes = [4,1], marker='+', markersize=3, color=oddColor, label=l)

    data = evenSpectrum["loc"][1]-2*evenSpectrum["loc"][0]
    l=label("loc", 1)
    plt.plot(xList, data, linewidth=0.7,
            dashes = [4,1], marker='+', markersize=3, color=evenColor, label=l)

#############################################################################

    plt.figure(1, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
    #plt.xlim(min(xList)-0.01, max(xList)+0.01)
    plt.title("m=1, dim0="+str(dim0)+", g="+str(g)+", n0eigs="+str(n0eigs)+", ordM="+str(ordMass))
    plt.xlabel(r"$L$")
    plt.ylabel(r"$\Lambda$")
    plt.legend(loc='lower right', prop={'size':10})
    plt.savefig("fig_vacVsL_g="+str(g)+"_n0eigs="+str(n0eigs)+"_dim0="+str(dim0)+"_ordM="+str(ordMass)+".svg")

    plt.figure(2, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
    plt.title("m=1, dim0="+str(dim0)+", g="+str(g)+", n0eigs="+str(n0eigs)+", ordM="+str(ordMass))
    plt.xlabel(r"$L$")
    plt.ylabel(r"$m$")
    plt.legend(loc='lower right', prop={'size':6})
    plt.savefig("fig_massVsL_g="+str(g)+"_n0eigs="+str(n0eigs)+"_dim0="+str(dim0)+"_ordM="+str(ordMass)+".svg")

    plt.figure(3, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
    plt.title("m=1, dim0="+str(dim0)+", g="+str(g)+", n0eigs="+str(n0eigs)+", ordM="+str(ordMass))
    plt.ylim(-0.03,-0.00)
    plt.xlabel(r"$L$")
    plt.ylabel(r"$E_2-2 E_1$")
    plt.legend(loc='lower right', prop={'size':6})
    plt.savefig("fig_boundStateVsL_g="+str(g)+"_n0eigs="+str(n0eigs)+"_dim0="+str(dim0)+"_ordM="+str(ordMass)+".svg")


if __name__ == "__main__":
    main(sys.argv)
