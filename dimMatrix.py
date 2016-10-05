import inspect
import os
import phi4
import sys
import math
import scipy
from statefuncs import *

k = 1
m = 1.

argv = sys.argv

args = "<L> <Emax> <Ebarmax>"
if len(argv) < 4:
    print("python", argv[0], args)
    sys.exit(-1)

L = float(argv[1])
Emax = float(argv[2])
Ebarmax = float(argv[3])

Ebarlist = scipy.linspace(Emax+1, Ebarmax, int(Ebarmax-Emax))

print(Ebarlist)

a = phi4.Phi4(m,L)
a.buildBasis(Emax=Emax)
a.buildMatrix(k)


nmax =  a.helper.nmax

subbasis = Basis(k, a.basis[k].stateList[:50], a.helper)


for Ebar in Ebarlist:

    a.computeVHH(k, subbasis, Emax, Ebar)

    # print("HE basis size", a.basisH[k].size)

    print("{",Ebar,",", len(a.VHH[k].M.data),",",a.nops,",",a.basisH[k].size,"},")
