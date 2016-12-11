import sys
import matplotlib.pyplot as plt
import scipy
import math
from scipy import pi, log, log10, array, sqrt, stats
from matplotlib import rc
from cycler import cycler
import database

ntails = {1:117, -1:108}

# List of all the contributions to DH3. Sequentially, we add DH3<<, DH3<> and DH3>>
# tlist = ((False,False,False),(True,True,True))
tlist = ((True,True,True),)

# Ratio between EL and ET
ratioELET = 3
# Ratio between ELpp and ELp
ratioELppELp = 1.5

neigs = 2

klist = (1,-1)

output = "png"
# renlist = ("raw", "renloc", "rentails")

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)



def plotvsELp(ELplist, tlist):

    xlist = ELplist

    db = database.Database()

    approxQuery = {"g":g, "L":L, "EL":EL, "ET":ET}

    spectrum = {k:{t:[] for t in tlist} for k in klist}
    masses = {k:{t:[] for t in tlist} for k in klist}

    for k in klist:

        for ELp in ELplist:
            approxQuery["ELp"] = ELp
            approxQuery["ELpp"] = ratioELppELp*ELp

            for t in tlist:
                nonloc3mix, loc3mix, loc3 = t
                exactQuery = {"loc3":loc3, "loc3mix":loc3mix, "nonloc3mix":nonloc3mix,
                    "ren":"rentails"}

                try:
                    exactQuery["k"] = 1
                    exactQuery["ntails"] = ntails[1]
                    spectrum[k][t].append(db.getObjList('spec', exactQuery, approxQuery)[0])

                    exactQuery["k"] = -1
                    exactQuery["ntails"] = ntails[-1]
                    spectrum[k][t].append(db.getObjList('spec', exactQuery, approxQuery)[0])

                except IndexError:
                    print("Not found: ", exactQuery, approxQuery)
                    sys.exit(-1)

    # Convert to array
    for k in klist:
        for t in tlist:
            spectrum[k][ren] = array(spectrum[k][ren])


    # Mass
    if -1 in klist and 1 in klist:
        for k in klist:
            for t in tlist:
                for i in range(neigs):
                    masses[k][t].append(spectrum[k][t][:,i]-spectrum[1][t][:,0])

    # Convert to array
    for k in klist:
        for t in tlist:
            masses[k][t] = array(masses[k][t])


    # SPECTRUM
    for k in klist:
        plt.figure(fignum(k))
        for i in range(neigs):
            for ren in renlist:
                data = spectrum[k][ren][:,i]
                if i==0:
                    label="ren="+ren
                else:
                    label = None
                plt.plot(xlist, data, label="ren="+ren)
            plt.gca().set_prop_cycle(None)

    # MASS
    if -1 in klist and 1 in klist:
        plt.figure(3)
        for k in (-1,1):
            for i in range(neigs):
                if k==1 and i==0:
                    continue
                for ren in renlist:
                    data = masses[k][ren][i]
                    if i==0:
                        label="ren="+ren
                    else:
                        label = None
                    plt.plot(xlist, data, label=label)
                plt.gca().set_prop_cycle(None)

argv = sys.argv


if len(argv) < 6:
    print(argv[0], "<L> <g> <ET> <ELpmin> <ELpmax>")
    sys.exit(-1)

L = float(argv[1])
g = float(argv[2])
ET = float(argv[3])
ELpmin = float(argv[4])
ELpmax = float(argv[5])


EL = ratioELET*ET

ELplist = scipy.linspace(ELpmin, ELpmax, (ELpmax-ELpmin)*2+1)
print("ELplist:", ELplist)


params = {'legend.fontsize': 8}
plt.rcParams.update(params)

plt.rc('axes', prop_cycle=(cycler('color', ['r', 'g', 'b', 'y']) +
    cycler('linestyle', ['-', '--', ':', '-.'])))


plotvsELp(ELplist, tlist)


title = r"$g$={0:.1f}, $L$={1:.1f}, $E_T$={2:.1f}, $E_L$={3:.1f},$E_L''={4} E_L'$".format(g,L,
        ET,EL,ratioELppELp)
fname = "g={0:.1f}_L={1:.1f}_ET={2:.1f}_EL={3:.1f}_ELpp/ELp={4}.{5}".format(g,L,ET,EL,
        ratioELppELp,output)
loc = "lower right"


# Even eigenvalues
if 1 in klist:
    plt.figure(1, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
    plt.title(title)
    plt.xlabel(r"$E_L'$")
    plt.ylabel(r"$E_i$ even")
    plt.legend(loc=loc)
    plt.savefig("evenvsELp_"+fname)


# Odd eigenvalues
if -1 in klist:
    plt.figure(2, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
    plt.title(title)
    plt.xlabel(r"$E_L'$")
    plt.ylabel(r"$E_i$ odd")
    plt.legend(loc=loc)
    plt.savefig("oddvsELp_"+fname)


# Mass
if 1 in klist and -1 in klist:
    plt.figure(3, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
    plt.title(title)
    plt.xlabel(r"$E_L'$")
    plt.ylabel(r"$m_{\rm ph}$")
    plt.legend(loc=loc)
    plt.savefig("massvsELp_"+fname)
