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

klist = (1,-1)

neigs = 3


# Ratio between EL and ET
ratioELET = 1.5
# Ratio between ELp and ET
ratioELpET = 1.5
# Ratio between ELpp and ELp
ratioELppELp = 1.5

maxntails = 200


def fignum(k):
    if k==1:
        return 1
    elif k==-1:
        return 2


def plotvsET(ETlist):

    xlist = ETlist

    db = database.Database()

    exactQuery = {"ren":"rentails"}
    approxQuery = {"g":g, "L":L}
    boundQuery = {"ntails": (0,maxntails)}

    spectrum = {k:[] for k in klist}
    ntails = {k:[] for k in klist}

    for k in klist:
        exactQuery["k"] = k

        for ET in ETlist:

            # print("ET=",ET)

            EL = ratioELET*ET
            ELp = ratioELpET*ET
            ELpp = ratioELppELp*ELp

            approxQuery["ET"] = ET
            approxQuery["EL"] = EL
            approxQuery["ELp"] = ELp
            approxQuery["ELpp"] = ELpp


            # NOTE We select the entry with the highest ntails below maxntails
            spectrum[k].append(db.getObjList('spec', exactQuery, approxQuery,
                boundQuery, orderBy="ntails")[-1])

            ntails[k].append(db.getObjList('ntails', exactQuery, approxQuery,
                boundQuery, orderBy="ntails")[-1])


    # SPECTRUM
    for k in klist:
        plt.figure(fignum(k))
        sp = array(spectrum[k])
        for i in range(neigs):
            data = sp[:,i]
            plt.plot(xlist, data)

    plt.figure(3)
    for k in klist:
        y = array(ntails[k])
        plt.plot(xlist,y,label="k="+str(k))


argv = sys.argv


if len(argv) < 5:
    print(argv[0], "<L> <g> <ETmin> <ETmax>")
    sys.exit(-1)

L = float(argv[1])
g = float(argv[2])
ETmin = float(argv[3])
ETmax = float(argv[4])

ETlist = scipy.linspace(ETmin, ETmax, (ETmax-ETmin)*2+1)
print("ETlist:", ETlist)


params = {'legend.fontsize': 8}
plt.rcParams.update(params)

plt.rc('axes', prop_cycle=(cycler('color', ['r', 'g', 'b', 'y']) +
    cycler('linestyle', ['-', '--', ':', '-.'])))


plotvsET(ETlist)



# NOTE plot of the spectrum

for k in klist:

    title = r"$g$={0:.1f}, $L$={1:.1f}, $k$={2}, maxntails={3},"\
                "$E_L/E_T$={4:.1f}, $E_L'/E_T$={5:.1f}, $E_L''/E_T={6:.1f}$"\
                .format(g,L,k,maxntails,ratioELET,ratioELpET,ratioELppELp)
    fname = "g={0:.1f}_L={1:.1f}_maxntails={2}.{3}_"\
                "ELET={3}_ELpET={4}_ELppELp={5}.{6}"\
                .format(g,L,maxntails,ratioELET,ratioELpET,ratioELppELp,output)
    loc = "lower right"

    plt.figure(fignum(k), figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
    plt.title(title)
    plt.xlabel(r"$E_{T}$")
    plt.ylabel(r"$E_i$")
# plt.legend(loc=loc)

    if k==1:
        plt.savefig("evenvsET_"+fname)
    else:
        plt.savefig("oddvsET_"+fname)




# NOTE Plot of the number of tails

title = r"$g$={0:.1f}, $L$={1:.1f}, maxntails={2}".format(g,L,maxntails)
fname = "g={0:.1f}_L={1:.1f}_maxntails={2}.{3}".format(g,L,maxntails,output)
loc = "lower right"

plt.figure(3, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
plt.title(title)
plt.xlabel(r"$E_{T}$")
plt.ylabel(r"numtails")
plt.legend(loc=loc)

plt.savefig("tailsvsET_"+fname)
