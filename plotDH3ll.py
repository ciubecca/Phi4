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
g = 1.
eps = {g: -1}

test = False

# (DH3<> nonloc, DH3<> loc, DH3>> loc)
tlist = ((False,False,False),(True,True,False),(True,True,True))

indexList = [(0,0), # V0
        (0,14) # V6 + V2V4
        ]

# Ratio between ELp and ET
ratioELpET = 2
# Ratio between ELpp and ELp
ratioELppELp = 1.5
frac = ratioELppELp

params = {'legend.fontsize': 8}
plt.rcParams.update(params)
plt.rc('axes', prop_cycle=(cycler('color', ['r', 'g', 'b', 'y']) +
    cycler('linestyle', ['-', '--', ':', '-.'])))
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
output = "pdf"

argv = sys.argv

args = "<L> <ET> <ELpmin> <ELpmax>"
if len(argv) < 5:
    print("python", argv[0], args)
    sys.exit(-1)

L = float(argv[1])
ET = float(argv[2])
ELpmin = float(argv[3])
ELpmax = float(argv[4])


ELplist = scipy.linspace(ELpmin, ELpmax, (ELpmax-ELpmin)+1)
print("ELplist:", ELplist)


a = phi4.Phi4(m,L,k)
a.buildBasis(Emax=ET)
a.computePotential()
a.setglist(glist=[g])

basisl = a.basis
# print(a.basis.stateList[0:20])
a.computeLEVs(basisl, loc3=True)

a.genHEBasis(EL=ELpmax, ELp=ELpmax, ELpp=ratioELppELp*ELpmax)
print("HE basis size", a.basisH.size)

print("Computing HE matrices")
a.computeHEVs()

ret = {t:{index:[] for index in indexList} for t in tlist}

for ELp in ELplist:
    ELpp = ratioELppELp *ELp
    print("ELp={}, ELpp={}".format(ELp,ELpp))

    a.calcVV3(ELp, eps, test=test)

    for t in tlist:
        nonloc3mix, loc3mix, loc3 = t

        DH3ll = a.computeDH3(ET=ET, ELp=ELp, ELpp=ELpp, eps=eps,
                loc3=loc3, loc3mix=loc3mix, nonloc3mix=nonloc3mix)[g].todense()

        for index in indexList:
            ret[t][index].append(DH3ll[index])


for i,index in enumerate(indexList):
    print("Index:{}, State:{}".format(index[1], a.basis.stateList[index[1]]))

    for t in tlist:
        if t == (False,False,False):
            label = "$\Delta H_3^{< <}$"
        elif t == (True,True,False):
            label = "$\Delta H_3^{< <} + \Delta H_3^{< >}$"
        elif t == (True,True,True):
            label = "$\Delta H_3^{< <} + \Delta H_3^{< >}+\Delta H_3^{> >}$"
        plt.plot(ELplist, ret[t][index], label=label)

    title = r"$g$={}, $L$={:.1f}, $E_T$={:.1f}, $E_L''={} E_L'$".format(g,L,ET,frac)
    fname = "DH3_g={}_L={:.1f}_ET={:.1f}_frac={}_index={}.{}".\
            format(g,L,ET,frac,str(index).replace(" ",""),output)
    loc = "lower right"

    state = ["0","6_0"]

    plt.figure(1, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
    plt.title(title)
    plt.xlabel(r"$E_{L}'$")
    plt.ylabel(r"$\langle {} | \Delta H_3 | {}\rangle$".\
            format(state[0],state[i]))
    plt.legend(loc=4)

    plt.savefig(fname)
    plt.clf()
