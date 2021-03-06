################
# Computation of vacuum energy in old fashioned perturbation theory at third order
#################


from sys import argv, exit
from math import factorial
import matplotlib.pyplot as plt
from statefuncs import *
import copy
from scipy import sparse
from time import time
import sys
from phi4 import *

m = 1

if len(argv) < 3:
    print("{} <L> <ELmax>".format(argv[0]))
    exit(1)


L = float(argv[1])
ELmax = float(argv[2])


print("L={}, ELmax={}".format(L, ELmax))

ET = 2
ELmin = 5
ELlist = np.linspace(ELmin, ELmax, 10)

subidx = {k:[0] for k in (1,-1)}
res = []
a = Phi4(m, L, ET)
basis = a.bases[1]

print("Generating high-energy basis...")
basisH = genHEBases(a.bases, subidx, ELmax, ELmax, V2=False, k=1)[1]

print("basis size = {}".format(basisH.size))

print("Generating low-high matrix...")
VlH = genVHl(basis, subidx[1], basisH, L)
print("VlH.shape: ", VlH.shape)

Vhh = genVhh(basisH, L)

checkSymmetric(Vhh)

prop = 1/(-np.array(basisH.energyList))


for EL in ELlist:
    idxlist = basisH.subidxlist(EL, Emin=2)
    VlHsub = subcolumns(VlH, idxlist)
    Vhhsub = submatrix(Vhh, idxlist)

    propsub = prop[idxlist]
    proj = scipy.sparse.spdiags(propsub, 0, len(propsub), len(propsub)).tocsc()

    deltaE = (VlHsub*proj*Vhhsub*proj*VlHsub.transpose()).todense()[0,0]

    res.append(deltaE)


vac = np.array(res)/(L**2)
plt.figure(1)
plt.plot(ELlist, vac)
plt.title("L={}".format(L))
plt.tight_layout()
plt.xlabel(r"$E_T$")
plt.ylabel(r"$\Delta_3 E_0$")
plt.savefig("vacpertCubic_L={}.pdf".format(L))
