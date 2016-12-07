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

neigs = 10

klist = (1,)

# Ratio between EL and ET
ratioELET = 1.5
# Ratio between ELp and ET
ratioELpET = 1.5
# Ratio between ELpp and ELp
ratioELppELp = 1.5


maxntails = 100
step = 10
startntails = 20
ntailsList = list(range(startntails, maxntails+step, step))
ntailsList = {1:ntailsList, -1:ntailsList}
print("ntailsList", ntailsList)


ETlist = [15,20,24]
print("ETlist", ETlist)


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
    spectrum = {}

    for k in klist:
        spectrum[k] = {ET:[] for ET in ETlist}

        xlist[k] = ntailsList[k]
        exactQuery["k"] = k

        for ET in ETlist:

            EL = ratioELET*ET
            ELp = ratioELpET*ET
            ELpp = ratioELppELp*ELp

            approxQuery = {"g":g, "L":L, "EL":EL, "ET":ET, "ELp":ELp, "ELpp":ELpp}


            for ntails in ntailsList[k]:
                exactQuery["ntails"] = ntails

                try:
                    spectrum[k][ET].append(db.getObjList('spec', exactQuery, approxQuery)[0])
                except IndexError as e:
                    print(approxQuery)
                    print(exactQuery)
                    raise(e)


    # SPECTRUM
    for k in klist:
        plt.figure(fignum(k))
        for ET in ETlist:
            data = array(spectrum[k][ET])[:,0]
            plt.plot(xlist[k], data, label="ET="+str(ET))



argv = sys.argv


if len(argv) < 3:
    print(argv[0], "<L> <g>")
    sys.exit(-1)

L = float(argv[1])
g = float(argv[2])


params = {'legend.fontsize': 8}
plt.rcParams.update(params)

plt.rc('axes', prop_cycle=(cycler('color', ['r', 'g', 'b', 'y']) +
    cycler('linestyle', ['-', '--', ':', '-.'])))


plotvsntails(ntailsList)


for k in klist:

    title = r"$g$={0:.1f}, $L$={1:.1f}, $E_L/E_T$={2:.1f}, $E_L'/E_T$={3:.1f},"\
            "$E_L''/E_L'$={4:.1f}".format(g,L,ratioELET,ratioELpET,ratioELppELp)
    fname = "g={0:.1f}_L={1:.1f}.{2}".format(g,L,output)
    loc = "upper right"

    plt.figure(fignum(k), figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
    plt.title(title)
    plt.xlabel("ntails")
    plt.ylabel(r"$E_i$")
    plt.legend(loc=loc)

    if k==1:
        plt.savefig("evenvsTailsAndET_"+fname)
    elif k==-1:
        plt.savefig("oddvsTailsAndET_"+fname)
