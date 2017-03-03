import sys
import matplotlib.pyplot as plt
import scipy
import math
from scipy import pi, log, log10, array, sqrt, stats
from matplotlib import rc
from cycler import cycler
import database
from sys import exit

output = "png"
renlist = ("raw", "rentails", "renloc")

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

klist = (1,-1)

neigs = 6

# Ratio between EL and ET
ratioELET = 3
# Ratio between ELp and ET
ratioELpET = 2
# Ratio between ELpp and ELp
ratioELppELp = 1.5


def fignum(k):
    if k==1:
        return 1
    elif k==-1:
        return 2

marker = 'o'
markersize = 2.5


def plotvsL(xlist):

    db = database.Database()

    spectrum = {k:{ren:[] for ren in renlist} for k in klist}
    masses = {k:{ren:[] for ren in renlist} for k in klist}
    ntails = {k:[] for k in klist}

    for k in klist:

        for L in xlist:

            for ren in renlist:

                EL = ratioELET*ET
                ELp = ratioELpET*ET
                ELpp = ratioELppELp*ELp

                approxQuery = {"g":g, "L":L, "ET":ET}
                exactQuery = {"k": k, "ren":ren, "neigs":neigs}
                boundQuery = {}

                if ren=="rentails":
                    exactQuery["maxntails"] = None
                    exactQuery["tailsComputedAtET"] = ET
                    approxQuery["EL"] = EL
                    approxQuery["ELp"] = ELp
                    approxQuery["ELpp"] = ELpp

                try:
                    spectrum[k][ren].append(db.getObjList('spec', exactQuery,
                        approxQuery, boundQuery, orderBy="date")[0])

                except IndexError:
                    print("Not found:", exactQuery, approxQuery)
                    exit(-1)



    # Convert to array
    for k in klist:
        for ren in renlist:
            spectrum[k][ren] = array(spectrum[k][ren])


    # Mass
    for k in (-1,1):
        for ren in renlist:
            for i in range(int((1+k)/2), neigs):
                masses[k][ren].append(spectrum[k][ren][:,i]-spectrum[1][ren][:,0])

    # Convert to array
    for k in klist:
        for ren in renlist:
            masses[k][ren] = array(masses[k][ren])

    # SPECTRUM
    for k in klist:
        plt.figure(fignum(k))
        # for i in range(neigs):
        for i in range(1):
            for ren in renlist:
                data = spectrum[k][ren][:,i]
                if i==0:
                    label="ren="+ren
                else:
                    label = None
                if ren == "renloc":
                    print("k=",k, " ", ",".join(str(x) for x in data))
                plt.plot(xlist, data, label=label, marker=marker, markersize=markersize)
            plt.gca().set_prop_cycle(None)

    # MASS
    plt.figure(3)
    # for k in (-1,1):
    for k in (-1,):
        # for i in range(neigs-int((1+k)/2)):
        for i in range(1):
            for ren in renlist:
                data = masses[k][ren][i]
                if i==0:
                    label="ren="+ren
                else:
                    label = None
                plt.plot(xlist, data, label=label,
                        markersize=markersize, marker=marker)
            plt.gca().set_prop_cycle(None)

argv = sys.argv


if len(argv) < 5:
    print(argv[0], "<ET> <g> <Lmin> <Lmax>")
    sys.exit(-1)

g = float(argv[1])
ET = float(argv[2])
Lmin = float(argv[3])
Lmax = float(argv[4])

Llist = scipy.linspace(Lmin, Lmax, (Lmax-Lmin)*2+1)
# print("ETlist:", ETlist)

print("g=", g)

params = {'legend.fontsize': 8}
plt.rcParams.update(params)

plt.rc('axes', prop_cycle=(cycler('color', ['r', 'g', 'b', 'y']) +
    cycler('linestyle', ['-', '--', ':', '-.'])))


plotvsL(Llist)


title = r"$g$={0:.1f}, $E_T$={1:.1f},"\
            "$E_L/E_T$={2:.1f}, $E_L'/E_T$={3:.1f}, $E_L''/E_L'={4:.1f}$"\
            .format(g,ET,ratioELET,ratioELpET,ratioELppELp)
fname = "g={0:.1f}_ET={1:.1f}_"\
            "ELET={2}_ELpET={3}_ELppELp={4}.{5}"\
            .format(g,ET,ratioELET,ratioELpET,ratioELppELp,output)
loc = "upper right"


# Even eigenvalues
if 1 in klist:
    plt.figure(1, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
    plt.title(title)
    plt.xlabel(r"$L$")
    plt.ylabel(r"$E_i$ even")
    plt.legend(loc=loc)
    plt.savefig("evenvsL_"+fname)


# Odd eigenvalues
if -1 in klist:
    plt.figure(2, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
    plt.title(title)
    plt.xlabel(r"$L$")
    plt.ylabel(r"$E_i$ odd")
    plt.legend(loc=loc)
    plt.savefig("oddvsL_"+fname)


# Mass
if 1 in klist and -1 in klist:
    plt.figure(3, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
    plt.title(title)
    plt.xlabel(r"$L$")
    plt.ylabel(r"$m_{\rm ph}$")
    plt.legend(loc=loc)
    plt.savefig("massvsL_"+fname)
