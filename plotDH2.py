import inspect
import os
import phi4
import sys
import math
import scipy
from statefuncs import *
from matplotlib import rc
import matplotlib.pyplot as plt
from cycler import cycler

k = +1
m = 1
g = 1
eps = -1

indexList = [(0,0),(0,1)]


argv = sys.argv

args = "<L> <ET> <ELmin> <ELmax>"
if len(argv) < 5:
    print("python", argv[0], args)
    sys.exit(-1)

L = float(argv[1])
ET = float(argv[2])
ELmin = float(argv[3])
ELmax = float(argv[4])


ELlist = scipy.linspace(ELmin, ELmax, (ELmax-ELmin)*2+1)
print("ELlist:", ELlist)

a = phi4.Phi4(m, L, k)
a.buildBasis(Emax=ET)
a.computePotential()
print("Full basis size: ", a.basis.size)
glist = [g]
a.setglist(glist=glist)

a.computeLEVs(a.basis, loc3=False)
a.genHEBasis(EL=ELmax, ELp=None, ELpp=None)
a.computeHEVs()

# Whether or not we include DH2>
tlist = (False, True)
ret = {}

for t in tlist:
    loc2 = t
    ret[t] = {index:[] for index in indexList}

    for EL in ELlist:
        print("EL={}".format(EL))

        DH2lL = a.computeDH2(ET, EL, eps={g:eps}, loc2=loc2)[0][g].todense()

        for index in indexList:
            ret[t][index].append(DH2lL[index])


for index in indexList:
    params = {'legend.fontsize': 8}
    plt.rcParams.update(params)

    plt.rc('axes', prop_cycle=(cycler('color', ['r', 'g', 'b', 'y']) +
        cycler('linestyle', ['-', '--', ':', '-.'])))

    rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
    rc('text', usetex=True)

    for t in tlist:
        label=r"$\Delta H_2^<$"
        if t==True:
            label += r"$ + \Delta H_2^>$"
        plt.plot(ELlist, ret[t][index], label=label)


    output = "pdf"
    title = r"$L$={:.1f}, $E_T$={:.1f}".format(L,ET,index)
    fname = "DH2_L={:.1f}_ET={:.1f}_index={}.{}".\
            format(L,ET,str(index).replace(" ",""),output)

    state = ["0", "2_0"]

    plt.figure(1, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
    plt.title(title)
    plt.xlabel(r"$E_{L}$")
    plt.ylabel(r"$\langle {} | \Delta H_2 | {} \rangle$"\
            .format(state[index[0]], state[index[1]]))
    plt.legend(loc=1)

    plt.savefig(fname, bbox_inches='tight')
    plt.clf()
