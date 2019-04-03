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
from nlo import *
from paramplots import *
m = 1

if len(argv) < 3:
    print("{} <L> <EL>".format(argv[0]))
    exit(1)


def ct0(ET, En=0):
    return -24**2*1/(96*(4*pi)**3)*((ET-En)-8*log((ET-En)/4)-16/(ET-En))



L = float(argv[1])
ELmax = int(argv[2])

ET = 2
ELmin = 10


ELlist = np.array(list(range(ELmin, ELmax+1)))

print("L={}, ELmax={}".format(L, ELmax))
print("ELlist: ", ELlist)

subidx = {k:[0] for k in (-1,1)}
E0 = {1:0, -1:m}

res = {k:[] for k in (-1,1)}

a = Phi4(m, L, ET)
bases = a.bases

basesH = genHEBases(bases, subidx, ELmax, ELmax)

for k in (-1,1):
    basisH = basesH[k]

    V = genVHl(bases[k], subidx[k], basisH, L)

    prop = 1/(E0[k]-np.array(basisH.energyList))

    for EL in ELlist:
        idxlist = basisH.subidxlist(EL, Emin=2)
        Vsub = subcolumns(V, idxlist)
        propsub = prop[idxlist]

        deltaE = np.einsum("ij,j,kj", Vsub.todense(), propsub, Vsub.todense())[0][0]

        res[k].append(deltaE)


print("Vacuum:", res[1])


vac = np.array(res[1])/(L**2)-ct0(ELlist)
plt.figure(1)
plt.plot(ELlist, vac)
plt.title("L={}".format(L))
plt.tight_layout()
plt.xlabel(r"$E_T$")
plt.ylabel(r"$\Delta E_0 - c_0(E_T)$")
plt.savefig("vacpert_L={}.pdf".format(L))


mass = np.array(res[-1])-np.array(res[1]) - ct2(ELlist)
massnlo = mass - (L**2)*(ct0(ELlist,1)- ct0(ELlist))
plt.figure(2)
plt.plot(ELlist, mass, label="loc")
plt.plot(ELlist, massnlo, label="nlo")
plt.tight_layout()
plt.title("L={}".format(L))
plt.legend()
plt.xlabel(r"$E_T$")
plt.ylabel(r"$\Delta (E_1-E_0) - c_2(E_T)$")
plt.savefig("masspert_L={}.pdf".format(L))
