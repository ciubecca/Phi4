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

indexList = [(0,0),(0,1),(1,1)]

maxntails = 2
minoverlap = 10**-2

frac = 3

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


a = phi4.Phi4(m,L)

a.buildBasis(Emax=ET)

a.computePotential(k)

a.setCouplings(0,0,g)
a.computeEigval(k, ET, "raw")

vectorlist = [state for i,state in enumerate(a.basis[k])
        if abs(a.eigenvectors["raw"][1][0][i]) > minoverlap][:maxntails]
print(sorted(occn(state) for state in vectorlist))
basisl = Basis(k, vectorlist, a.basis[k].helper)
print("subbasis size:", basisl.size)

print(basisl)

a.genHEBases(k, basisl, EL=ELmax, ELpp=None)
print("HE basis size", a.basisH[k].size)

a.computeLEVs(k)

print("Computing HE matrices")
a.computeHEVs(k)

eps = -1

tlist = (False, True)

ret = {}

for t in tlist:
    loc2 = t
    ret[t] = {index:[] for index in indexList}

    for EL in ELlist:
        print("EL={}".format(EL))

        DH2lL = a.computeDH2(basisl, k, ET, EL, eps=eps, loc2=loc2).todense()

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
        plt.plot(ELlist, ret[t][index], label=str(t))

    output = "png"
    title = r"$L$={:.1f}, $E_T$={:.1f}, index={}".format(L,ET,index)
    fname = "DH2_L={:.1f}_ET={:.1f}_index={}.{}".format(L,ET,str(index),output)
    loc = "lower right"

    plt.figure(1, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
    plt.title(title)
    plt.xlabel(r"$E_{L}$")
    plt.ylabel(r"$\Delta H_2 {}$".format(index))
    plt.legend(loc=loc)


    plt.savefig(fname)

    plt.clf()
