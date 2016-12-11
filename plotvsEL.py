import sys
import matplotlib.pyplot as plt
import scipy
import math
from scipy import pi, log, log10, array, sqrt, stats
from matplotlib import rc
from cycler import cycler
import database

output = "png"
# renlist = ("raw", "renloc", "rentails")
renlist = ("rentails",)

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

neigs = 2

minoverlap = 10**(-2)
# Ratio between ELpp and ELp
ratio3 = 1.5


def plotvsEL(ELlist):

    xlist = ELlist
    nonloc3mix, loc3mix, loc3 = (True, True, True)

    db = database.Database()

    exactQuery = {"loc3":loc3, "loc3mix":loc3mix, "nonloc3mix":nonloc3mix,
            "ren":"rentails", "minoverlap":minoverlap}
    approxQuery = {"g":g, "L":L, "ELp":ELp, "ELpp":ELpp, "ET":ET}

    oddSp = []
    evenSp = []

    for EL in ELlist:
        approxQuery["EL"] = EL

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


if len(argv) < 7:
    print(argv[0], "<L> <g> <ET> <ELp> <ELmin> <ELmax>")
    sys.exit(-1)

L = float(argv[1])
g = float(argv[2])
ET = float(argv[3])
ELp = float(argv[4])
ELmin = float(argv[5])
ELmax = float(argv[6])

ELpp = ELp*ratio3

ELlist = scipy.linspace(ELmin, ELmax, (ELmax-ELmin)*2+1)
print("ELlist:", ELlist)


params = {'legend.fontsize': 8}
plt.rcParams.update(params)

plt.rc('axes', prop_cycle=(cycler('color', ['r', 'g', 'b', 'y']) +
    cycler('linestyle', ['-', '--', ':', '-.'])))


plotvsEL(ELlist)


title = r"$g$={0:.1f}, $L$={1:.1f}, $E_T$={2:.1f}, $E_L'$={3:.1f},$E_L''$={4:.1f}".format(g,L,ET,ELp,ELpp)
fname = "g={0:.1f}_L={1:.1f}_ET={2:.1f}_ELp={3:.1f}_ELpp={4:.1f}.{5}".format(g,L,ET,ELp,ELpp,output)
loc = "lower right"

plt.figure(1, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
plt.title(title)
plt.xlabel(r"$E_{L}$")
plt.ylabel(r"$E_i$ even")
# plt.legend(loc=loc)


plt.savefig("evenSp_"+fname)



plt.figure(2, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
plt.title(title)
plt.xlabel(r"$E_{L}$")
plt.ylabel(r"$E_i$ odd")
# plt.legend(loc=loc)


plt.savefig("oddSp_"+fname)
