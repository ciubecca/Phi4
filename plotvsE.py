import sys
import matplotlib.pyplot as plt
import scipy
import math
from scipy import pi, log, log10, array, sqrt, stats
from matplotlib import rc
import database

output = "pdf"
renlist = ("raw", "ren")

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

k = 1
neigs = 10
ETmin = 10
ETmax = 15
ELETdiff = 5

def plotvsE(Elist, ELETdiff):

    db = database.Database()

    exactQuery = {"k":k}
    approxQuery = {"g":g, "L":L}

    E0 = {}
    for ren in renlist:
        E0[ren] = []

        for ET in ETlist:
            exactQuery["ren"] = ren
            approxQuery["ET"] = ET
            if ren=="ren":
                approxQuery["EL"] = ET+ELETdiff
            E0[ren].append(db.getObjList('spec', exactQuery, approxQuery)[0][0])


    # VACUUM ENERGY
    plt.figure(1)

    linewidth = 1
    dashes = [4,4]
    marker = None
    markersize = 1

    data = E0["raw"]
    plt.plot(Elist, data, linewidth=linewidth, color="b", marker=marker,
            markersize=markersize, dashes = dashes, label="raw")

    # if tails:
        # plt.axhline(y=data[-1], color='k')
        # plt.axhline(y=data[-1], xmin=min(Elist), xmax=max(Elist), linewidth=2, color = 'k')

    data = E0["ren"]
    plt.plot(Elist, data, linewidth=linewidth, color="r", marker=marker,
            markersize=markersize, dashes = dashes, label="ren")



argv = sys.argv
if len(argv) < 6:
    print(argv[0], "<L> <g> <ETmin> <ETmax> <EL-ET>")
    sys.exit(-1)

L = float(argv[1])
g = float(argv[2])
ETmin = float(argv[3])
ETmax = float(argv[4])
ELETdiff = float(argv[5])

ETlist = scipy.linspace(ETmin, ETmax, ETmax-ETmin+1)
print("ETlist:", ETlist)


params = {'legend.fontsize': 8}
plt.rcParams.update(params)

plotvsE(ETlist, ELETdiff)

plt.figure(1, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
#plt.xlim(min(xList)-0.01, max(xList)+0.01)
plt.title(r"$g$={0:.1f}, $L$={1:.1f}, $E_L-E_T$={2:.1f}".format(g,L,ELETdiff))
plt.xlabel(r"$E_{{\rm max}}$")
plt.ylabel(r"$E_0$")
plt.legend(loc="lower right")

plt.savefig("figs/fig_E0vsET_g={0:.1f}_L={1:.1f}_ELETdiff={2:.1f}.{3}"
        .format(g,L,ELETdiff,output))



