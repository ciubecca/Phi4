from nlo import *
from sys import argv, exit
from math import factorial
import matplotlib.pyplot as plt
from statefuncs import *
import copy
from scipy import sparse
from database import *
from time import time
from phi4 import *
from paramplots import *
m = 1

if len(argv) < 4:
    print("{} <L> <EL> <Lambda>".format(argv[0]))
    exit(1)


def ct0(ET):
    return -24**2*1/(96*(4*pi)**3)*(ET-8*log(ET))

def ct2(ET):
    return -2/((4*pi)**2)*log(ET)

L = float(argv[1])
ELmax = float(argv[2])
Lambda = float(argv[3])

ET = 2
ELmin = 5


ELlist = np.linspace(ELmin, ELmax, 10)

print("L={}, Lambda={}, ELmax={}".format(L, Lambda, ELmax))
print("ELlist: ", ELlist)

bases = Basis.fromScratch(m, L, ET, Lambda)


subidx = [0]

E0 = {1:0, -1:m}

res = {k:[] for k in (-1,1)}

# for k in (-1,1):
for k in (1,):
    print("Computing basis...")
    basisH = genHEBasis(bases[k], subidx, ELmax, ELmax)
    print("k={} basis size={}".format(k, basisH.size))
    print("Computing matrix...")
    V = genVHl(bases[k], subidx, basisH, L)
    print("Done")
    prop = 1/(E0[k]-array(basisH.energyList))

    for EL in ELlist:
        idxlist = basisH.subidxlist(EL, Lambda)
        Vsub = subcolumns(V, idxlist)
        propsub = prop[idxlist]

        deltaE = np.einsum("ij,j,kj", Vsub.todense(), propsub, Vsub.todense())[0][0]
        res[k].append(deltaE)

vac = array(res[1])/(L**2)-ct0(ELlist)
# vac = array(res[1])/(L**2)
plt.plot(ELlist, vac)
plt.title("L={}".format(L))
plt.xlabel(r"$E_T$")
plt.ylabel(r"$\Delta E_0 - c_0(E_T)$")
plt.savefig("vacpert_L={}.pdf".format(L))

print(vac)
