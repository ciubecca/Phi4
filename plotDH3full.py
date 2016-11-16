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
eps = -1

indexList = [(0,0),(0,1),(1,1)]

# Ratio between ELpp and ELp
frac = 3

argv = sys.argv

args = "<L> <ET> <ELpmin> <ELpmax>"
if len(argv) < 5:
    print("python", argv[0], args)
    sys.exit(-1)

L = float(argv[1])
ET = float(argv[2])
ELpmin = float(argv[3])
ELpmax = float(argv[4])


ELplist = scipy.linspace(ELpmin, ELpmax, (ELpmax-ELpmin)*2+1)
print("ELplist:", ELplist)


a = phi4.Phi4(m,L)

a.buildBasis(Emax=ELpmax)

a.computePotential(k)

tlist = ((False,False,False),)

if True in (t[2] for t in tlist):
    a.calcVV3(ELplist, eps)

ret = {}

basis = a.basis[k]
energyArr = array(basis.energyList)
propagator = scipy.sparse.spdiags(1/(eps-energyArr), 0, basis.size, basis.size)
Vfull = a.V[k][4].M

for t in tlist:
    nonloc3mix, loc3mix, loc3 = t
    ret[t] = {index:[] for index in indexList}

    for ELp in ELplist:
        print("ELp={}".format(ELp))

        projh = scipy.sparse.spdiags(array([int(ET<e<ELp) for e in energyArr]), 0,
            basis.size,basis.size)

        DH3Full = Vfull*propagator*projh*Vfull*propagator*projh*Vfull

        for index in indexList:
            ret[t][index].append(DH3Full[index])


for index in indexList:
    params = {'legend.fontsize': 8}
    plt.rcParams.update(params)

    plt.rc('axes', prop_cycle=(cycler('color', ['r', 'g', 'b', 'y']) +
        cycler('linestyle', ['-', '--', ':', '-.'])))

    rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
    rc('text', usetex=True)

    for t in tlist:
        plt.plot(ELplist, ret[t][index], label=str(t))

    output = "png"
    title = r"$L$={:.1f}, $E_T$={:.1f}, $E_L''={} E_L'$, index={}".format(L,ET,frac,index)
    fname = "DH3Full_L={:.1f}_ET={:.1f}_frac={}_index={}.{}".format(L,ET,frac,str(index),output)
    loc = "lower right"

    plt.figure(1, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
    plt.title(title)
    plt.xlabel(r"$E_{L}'$")
    plt.ylabel(r"$\Delta H_3 {}$".format(index))
    plt.legend(loc=loc)


    plt.savefig(fname)

    plt.clf()
