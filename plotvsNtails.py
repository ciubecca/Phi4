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

klist = (1,-1)

# Ratio ELpp/ELp
ratio3 = 1.5

maxntails = 100
step = 10
startntails = 20
ntailsList = list(range(startntails, maxntails+step, step))
ntailsList = {1:ntailsList, -1:ntailsList}

print("ntailsList", ntailsList)

def fignum(k):
    if k==1:
        return 1
    elif k==-1:
        return 2


def plotvsntails(ntailsList):

    xlist = {}
    nonloc3mix, loc3mix, loc3 = True, True, True

    db = database.Database()

    exactQuery = {"loc3":loc3, "loc3mix":loc3mix, "nonloc3mix":nonloc3mix,
            "ren":"rentails"}
    approxQuery = {"g":g, "L":L, "EL":EL, "ET":ET, "ELp":ELp}
    approxQuery["ELpp"] = ELpp
    spectrum = {k:[] for k in klist}

    for k in klist:
        xlist[k] = ntailsList[k]
        exactQuery["k"] = k

        for ntails in ntailsList[k]:
            exactQuery["ntails"] = ntails

            spectrum[k].append(db.getObjList('spec', exactQuery, approxQuery)[0])




    # SPECTRUM
    for k in klist:
        plt.figure(fignum(k))
        sp = array(spectrum[k])
        for i in range(neigs):
            data = sp[:,i]
            plt.plot(xlist[k], data)



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


for k in klist:

    title = r"$g$={0:.1f}, $L$={1:.1f}, $E_T$={2:.1f}, $E_L$={3:.1f},$E_L'$={4:.1f},$E_L''$={5:.1f}".format(g,L,ET,EL,ELp,ELpp)
    fname = "g={0:.1f}_L={1:.1f}_ET={2:.1f}_EL={3:.1f}_ELp={4:.1f}_ELpp={5:.1f}.{6}".format(g,L,ET,EL,ELp,ELpp,output)
    loc = "lower right"

    plt.figure(fignum(k), figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
    plt.title(title)
    plt.xlabel("ntails")
    plt.ylabel(r"$E_i$")
    # plt.legend(loc=loc)

    if k==1:
        plt.savefig("evenvsTails_"+fname)
    elif k==-1:
        plt.savefig("oddvsTails_"+fname)
