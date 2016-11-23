import sys
import matplotlib.pyplot as plt
import scipy
import math
from scipy import pi, log, log10, array, sqrt, stats
from matplotlib import rc
from cycler import cycler
import database


output = "png"
renlist = ("raw", "renloc", "rentails")

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

neigs = 1

# Ratio ELpp/ELp
ratio3 = 1.5
ntailsList = range(2, 46, 4)

print("ntailsList", ntailsList)



def plotvsntails(ntailsList):

    xlist = ntailsList
    nonloc3mix, loc3mix, loc3 = True, True, True

    db = database.Database()

    exactQuery = {"loc3":loc3, "loc3mix":loc3mix, "nonloc3mix":nonloc3mix,
            "ren":"rentails"}
    approxQuery = {"g":g, "L":L, "EL":EL, "ET":ET, "ELp":ELp}
    approxQuery["ELpp"] = ELpp
    oddSp = []
    evenSp = []

    for ntails in ntailsList:
        exactQuery["ntails"] = ntails

        exactQuery["k"] = 1
        evenSp.append(db.getObjList('spec', exactQuery, approxQuery)[0])

        exactQuery["k"] = -1
        oddSp.append(db.getObjList('spec', exactQuery, approxQuery)[0])



    evenSp = array(evenSp)
    oddSp = array(oddSp)

    # EVEN SPECTRUM
    plt.figure(1)

    for i in range(neigs):
        data = evenSp[:,i]
        plt.plot(xlist, data)


    # ODD SPECTRUM
    plt.figure(2)

    for i in range(neigs):
        data = oddSp[:,i]
        plt.plot(xlist, data)


argv = sys.argv


if len(argv) < 6:
    print(argv[0], "<L> <g> <ET> <EL> <ELp>")
    sys.exit(-1)

L = float(argv[1])
g = float(argv[2])
ET = float(argv[3])
EL = float(argv[4])
ELp = float(argv[5])

ELpp = ratio3*ELp


params = {'legend.fontsize': 8}
plt.rcParams.update(params)

plt.rc('axes', prop_cycle=(cycler('color', ['r', 'g', 'b', 'y']) +
    cycler('linestyle', ['-', '--', ':', '-.'])))


plotvsntails(ntailsList)


title = r"$g$={0:.1f}, $L$={1:.1f}, $E_T$={2:.1f}, $E_L$={3:.1f},$E_L'$={4:.1f},$E_L''$={5:.1f}$".format(g,L,ET,EL,ELp,ELpp)
fname = "g={0:.1f}_L={1:.1f}_ET={2:.1f}_EL={3:.1f}_ELp={4:.1f}_ELpp={5:.1f}.{6}".format(g,L,ET,EL,ELp,ELpp,output)
loc = "lower right"

plt.figure(1, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
plt.title(title)
plt.xlabel("ntails")
plt.ylabel(r"$E_i$ even")
plt.legend(loc=loc)


plt.savefig("figs/evenSp_"+fname)



plt.figure(2, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
plt.title(title)
plt.xlabel("ntails")
plt.ylabel(r"$E_i$ odd")
plt.legend(loc=loc)


plt.savefig("figs/oddSp_"+fname)
