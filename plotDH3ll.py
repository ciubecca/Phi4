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
g = 2.
eps = {g: -1}

test = False

# (DH3<> nonloc, DH3<> loc, DH3>> loc)
tlist = ((True,True,False),(True,True,True))

vectorlist = [[],
        [(0,2)],
        [(0,4)],
        [(0,6)],
        [(0,8)],
        [(0,1),(1,4),(-4,1)]
        ]
indexList = [(0,0), # V0
        (0,1), # V2
        (0,2), # V4 + V0V4
        (0,3), # V6 + V2V4
        (0,4), # V4V4
        (0,5) # V6
        ]


# Ratio between ELpp and ELp
ratioELppELp = 1.5
frac = ratioELppELp

params = {'legend.fontsize': 8}
plt.rcParams.update(params)
plt.rc('axes', prop_cycle=(cycler('color', ['r', 'g', 'b', 'y']) +
    cycler('linestyle', ['-', '--', ':', '-.'])))
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
output = "png"

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


a = phi4.Phi4(m,L)
a.buildBasis(Emax=ET)
a.computePotential(k)
a.setglist(glist=[g])


basisl = a.basis[k]
a.computeLEVs(k, basisl, loc3=True)

a.genHEBasis(k, EL=ELpmax, ELp=ELpmax, ELpp=ratioELppELp*ELpmax)
print("HE basis size", a.basisH[k].size)


print("Computing HE matrices")
a.computeHEVs(k)

ret = {t:{index:[] for index in indexList} for t in tlist}

for ELp in ELplist:
    ELpp = ratioELppELp *ELp
    print("ELp={}, ELpp={}".format(ELp,ELpp))

    a.calcVV3(ELp, eps, test=test)

    for t in tlist:
        nonloc3mix, loc3mix, loc3 = t

        DH3ll = a.computeDH3(k, ET=ET, ELp=ELp, ELpp=ELpp, eps=eps,
                loc3=loc3, loc3mix=loc3mix, nonloc3mix=nonloc3mix)[g].todense()

        for index in indexList:
            ret[t][index].append(DH3ll[index])


for index in indexList:

    for t in tlist:
        plt.plot(ELplist, ret[t][index], label=str(t))

    title = r"$g$={}, $L$={:.1f}, $E_T$={:.1f}, $E_L''={} E_L'$, index={}".\
            format(g,L,ET,frac,index)
    fname = "DHplots/DH3_g={}_L={:.1f}_ET={:.1f}_frac={}_index={}.{}".\
            format(g,L,ET,frac,str(index),output)
    loc = "lower right"

    plt.figure(1, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
    plt.title(title)
    plt.xlabel(r"$E_{L}'$")
    plt.ylabel(r"$\Delta H_3 {}$".format(index))
    plt.legend(loc=loc)

    plt.savefig(fname)

    plt.clf()
