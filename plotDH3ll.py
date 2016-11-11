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

minoverlap = 10**-2

def ELppf(ELp):
    return 1.5*ELp


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

a.buildBasis(Emax=ET)

a.computePotential(k)

a.setCouplings(0,0,g)
a.computeEigval(k, ET, "raw")

vectorlist = [state for i,state in enumerate(a.basis[k])
        if abs(a.eigenvectors["raw"][1][0][i]) > minoverlap]
print(sorted(occn(state) for state in vectorlist))
basisl = Basis(k, vectorlist, a.basis[k].helper)
print("subbasis size:", basisl.size)


a.genHEBases(k, basisl, EL=None, ELpp=ELppf(ELpmax))
print("HE basis size", a.basish[k].size)

a.computeLEVs(k)

print("Computing HE matrices")
a.computeHEVs(k)

eps = -1

ret = {}
index = (0,0)
tlist = ((False,False,False),(True,False,False),(True,True,False))

for t in tlist:
    nonloc3mix, loc3mix, loc3 = t
    ret[t] = []

    for ELp in ELplist:
        ELpp = ELppf(ELp)
        print("ELp={}, ELpp={}".format(ELp,ELpp))

        DH3ll = a.computeDH3(basisl, k, ET, ELp, ELpp=ELpp, eps=eps,
                loc3=loc3, loc3mix=loc3mix, nonloc3mix=nonloc3mix).todense()

        ret[t].append(DH3ll[index])

params = {'legend.fontsize': 8}
plt.rcParams.update(params)

plt.rc('axes', prop_cycle=(cycler('color', ['r', 'g', 'b', 'y']) +
    cycler('linestyle', ['-', '--', ':', '-.'])))

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

for t in tlist:
    plt.plot(ELplist, ret[t], label=str(t))


output = "png"
title = r"$L$={:.1f}, $E_T$={:.1f}, $E_L''=1.5 E_L'$, index={}".format(L,ET,index)
fname = "DH3_L={:.1f}_ET={:.1f}_index={}.{}".format(L,ET,str(index),output)
loc = "lower right"

plt.figure(1, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
plt.title(title)
plt.xlabel(r"$E_{L}'$")
plt.ylabel(r"$\Delta H_3$")
plt.legend(loc=loc)


plt.savefig(fname)
