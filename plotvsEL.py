import sys
import matplotlib.pyplot as plt
import scipy
import math
from scipy import pi, log, log10, array, sqrt, stats
from matplotlib import rc
from cycler import cycler
import database
from sys import exit
import numpy as np

dbname = "data/spectravsEL.db"
test = False

output = "pdf"

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

params = {'legend.fontsize': 8, 'lines.markersize':2.5, 'lines.marker':"o"}
plt.rcParams.update(params)

plt.style.use('ggplot')

plt.rc('axes', prop_cycle=(
    cycler('marker', ['x', 'o', 'v'])
    +cycler('linestyle', ['-', '--', ':'])
    +cycler('markersize', [5.,3.,3.])
    ))

color = {1:"b", -1:"r"}

k = 1
neigs = 6

# Ratio between ELp and ET
ratioELpET = 2
# Ratio between ELpp and ELp
ratioELppELp = 1.5


def plotvsEL(ELlist):

    db = database.Database(dbname)

    for loc2 in [False,True]:
        spectrum = np.array([db.getEigs(k,"rentails",g,L,ET,EL=EL,loc2=loc2,
            test=test) for EL in ELlist]).transpose()

        # VACUUM
        plt.figure(1)

        label=r"$\Delta H_2^<$"
        if loc2==True:
            label += r"$ + \Delta H_2^>$"

        data = spectrum[0]
        plt.plot(ELlist, data, label=label, color='k')


    # SPECTRUM
#     plt.figure(2)
    # for k, n0 in ((1,1),(-1,0)):

        # for i in range(n0,neigs):

            # for ren in renlist:
                # data = spectrum[k][ren][i]-spectrum[1][ren][0]
                # if i==n0:
                    # label="k={}, ren={}".format(k,ren)
                # else:
                    # label = None
                # plt.plot(ETlist[ren], data, label=label, color=color[k])

            # plt.gca().set_prop_cycle(None)

    # MASS
    # plt.figure(3)
    # for ren in renlist:
        # data = spectrum[-1][ren][0]-spectrum[1][ren][0]
        # label="ren={}".format(ren)
        # plt.plot(ETlist[ren], data, label=label, color="k")
        # ymin[-1] = min(ymin[-1], min(data))
        # ymax[-1] = max(ymax[-1], max(data))

    plt.gca().set_prop_cycle(None)


argv = sys.argv


if len(argv) < 6:
    print(argv[0], "<L> <g> <ET> <ELmin> <ELmax>")
    sys.exit(-1)

L = float(argv[1])
g = float(argv[2])
ET = float(argv[3])
ELmin = float(argv[4])
ELmax = float(argv[5])


ELlist = scipy.linspace(ELmin, ELmax, (ELmax-ELmin)*2+1)
print("ELlist:", ELlist)

plotvsEL(ELlist)


title = r"$g$={0:.1f}, $L$={1:.1f}".format(g,L)
fname = "g={0:.1f}_L={1:.1f}.{2}".format(g,L,output)

# VACUUM
plt.figure(1, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
plt.title(title)
plt.xlabel(r"$E_{L}$")
plt.ylabel(r"$\mathcal{E}_0$")
dx = 5*10**-1
plt.legend(loc=1)
plt.savefig("E0vsEL_"+fname, bbox_inches='tight')


# MASSES
# plt.figure(2, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
# plt.title(title)
# plt.xlabel(r"$E_{T}$")
# plt.ylabel(r"$\mathcal{E}_I - \mathcal{E}_0$")
# plt.xlim(xlims)
# plt.legend(loc=2)
# plt.savefig("specvsET_"+fname)

# MASS
# plt.figure(3, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
# plt.title(title)
# plt.xlabel(r"$E_{T}$")
# plt.ylabel(r"$\mathcal{E}_1 - \mathcal{E}_0$")
# plt.xlim(xlims[0]-dx, xlims[1]+dx)
# plt.ylim(ymin[-1]-10**-3,ymax[-1]+10**-3)
# plt.legend(loc=1)
# plt.savefig("massvsET_"+fname, bbox_inches='tight')
