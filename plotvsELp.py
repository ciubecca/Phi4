import sys
import matplotlib.pyplot as plt
import scipy
import math
from scipy import pi, log, log10, array, sqrt, stats
from matplotlib import rc
from cycler import cycler
import database


minoverlap = 3*10**(-3)
# List of all the contributions to DH3. Sequentially, we add DH3<<, DH3<> and DH3>>
tlist = ((False,False,False),(True,True,False),(True,True,True))
# Ratio between ELpp and ELp
ratio = 2


output = "png"
# renlist = ("raw", "renloc", "rentails")

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

neigs = 3


def ELppf(ELp):
    return 1.5*ELp


def plotvsELp(ELplist, t):

    xlist = ELplist
    nonloc3mix, loc3mix, loc3 = t

    db = database.Database()

    exactQuery = {"loc3":loc3, "loc3mix":loc3mix, "nonloc3mix":nonloc3mix,
            "ren":"rentails"}
    approxQuery = {"g":g, "L":L, "EL":EL, "ET":ET, "minoverlap":minoverlap}

    oddSp = []
    evenSp = []

    for ELp in ELplist:
        approxQuery["ELp"] = ELp
        approxQuery["ELpp"] = ratio*ELp

        exactQuery["k"] = 1
        evenSp.append(db.getObjList('spec', exactQuery, approxQuery)[0])

        exactQuery["k"] = -1
        oddSp.append(db.getObjList('spec', exactQuery, approxQuery)[0])


    nonloc3mix, loc3mix, loc3 = t
    label = str(t)

    evenSp = array(evenSp)
    oddSp = array(oddSp)

    # EVEN SPECTRUM
    plt.figure(1)

    for i in range(neigs):
        data = evenSp[:,i]
        plt.plot(xlist, data, label=label)


    # ODD SPECTRUM
    plt.figure(2)

    for i in range(neigs):
        data = oddSp[:,i]
        plt.plot(xlist, data, label=label)


argv = sys.argv


if len(argv) < 6:
    print(argv[0], "<L> <g> <ET> <ELpmin> <ELpmax> [<EL>]")
    sys.exit(-1)

L = float(argv[1])
g = float(argv[2])
ET = float(argv[3])
ELpmin = float(argv[4])
ELpmax = float(argv[5])

try:
    EL = float(argv[6])
except IndexError:
    EL = ratio*ELpmax
# ELETdiff = float(argv[5])


ELplist = scipy.linspace(ELpmin, ELpmax, (ELpmax-ELpmin)*2+1)
print("ELplist:", ELplist)


params = {'legend.fontsize': 8}
plt.rcParams.update(params)

plt.rc('axes', prop_cycle=(cycler('color', ['r', 'g', 'b', 'y']) +
    cycler('linestyle', ['-', '--', ':', '-.'])))



for t in tlist:
    plotvsELp(ELplist, t)


title = r"$g$={0:.1f}, $L$={1:.1f}, $E_T$={2:.1f}, $E_L$={3:.1f},$E_L''={4} E_L'$".format(g,L,ET,EL,ratio)
fname = "g={0:.1f}_L={1:.1f}_ET={2:.1f}_EL={3:.1f}_ratio={4}.{5}".format(g,L,ET,EL,ratio,output)
loc = "lower right"

plt.figure(1, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
plt.title(title)
plt.xlabel(r"$E_{L}'$")
plt.ylabel(r"$E_i$ even")
plt.legend(loc=loc)


plt.savefig("figs/evenSp_"+fname)



plt.figure(2, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
plt.title(title)
plt.xlabel(r"$E_{L}'$")
plt.ylabel(r"$E_i$ odd")
plt.legend(loc=loc)


plt.savefig("figs/oddSp_"+fname)
