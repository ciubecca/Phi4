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

if len(argv) < 3:
    print("{} <L> <Lambdamax>".format(argv[0]))
    exit(1)


# def ct0(Lambda, En=0):
    # coeff =
    # return -24**2*1/(96*(4*pi)**3)*((ET-En)-8*log((ET-En)/4)-16/(ET-En))


def ct2(Lambda):
    a = 1/(12*(4*pi)**2)
    b = 3.736124473715983
    c = -a * 1.5848415795962967
    return (24**2)*(-a*log(Lambda/(b*m)) + c*m/Lambda)


L = float(argv[1])
Lambdamax  = float(argv[2])

Lambdamin = 4
# This is enough to reproduce the full mass and vacuum corrections
ET = 4*Lambdamax+6

lamlist = np.linspace(Lambdamin, Lambdamax, 10)

print("L={}, Lambdamax={}, ET={}".format(L, Lambdamax, ET))

bases = Basis.fromScratch(m, L, Emax=2, Lambda=Lambdamax)

subidx = {k:[0] for k in (-1,1)}

E0 = {1:0, -1:m}

res = {k:[] for k in (-1,1)}

helper = None

print("Computing basis...")
basesH = genHEBases(bases, subidx, ET, ET)

for k in (-1,1):
    # basisH = genHEBasis(bases[k], subidx, ET, ET, helper)
    # helper = basisH.helper

    basisH = basesH[k]

    print("k={} basis size={}".format(k, basisH.size))
    print("Computing matrix...")
    V = genVHl(bases[k], subidx, basisH, L)
    print("Done")
    prop = 1/(E0[k]-array(basisH.energyList))

    for lam in lamlist:
        idxlist = basisH.subidxlist(ET, lam)
        Vsub = subcolumns(V, idxlist)
        propsub = prop[idxlist]

        deltaE = np.einsum("ij,j,kj", Vsub.todense(), propsub, Vsub.todense())[0][0]
        res[k].append(deltaE)



# vac = array(res[1])/(L**2)-ct0(ELlist)
# plt.figure(1)
# plt.plot(ELlist, vac)
# plt.title("L={}".format(L))
# plt.tight_layout()
# plt.xlabel(r"$E_T$")
# plt.ylabel(r"$\Delta E_0 - c_0(E_T)$")
# plt.savefig("vacpert_L={}.pdf".format(L))

mass = array(res[-1])-array(res[1]) - ct2(lamlist)
plt.figure(2)
plt.plot(lamlist, mass)
plt.tight_layout()
plt.title("L={}, ET={}".format(L, ET))
plt.xlabel(r"$\Lambda$")
plt.ylabel(r"$\Delta (E_1-E_0) - c_2(\Lambda)$")
plt.savefig("masspertvsLam_L={}.pdf".format(L))


# print(vac)
print(mass)

