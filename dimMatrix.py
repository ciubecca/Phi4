import inspect
import os
import phi4
import sys
import math
import scipy
import time
from statefuncs import *

k = 1
m = 1.

argv = sys.argv

args = "<L> <ET> <ELmax>"
if len(argv) < 4:
    print("python", argv[0], args)
    sys.exit(-1)

L = float(argv[1])
ET = float(argv[2])
ELmax = float(argv[3])

ELlist = scipy.linspace(ET+1, ELmax, int(ELmax-ET))

print(ELlist)

a = phi4.Phi4(m,L)
a.buildBasis(Emax=ET)

subbasis = Basis(k, a.basis[k].stateList[:50], a.basis[k].helper)

print("{Ebar, len(a.VHH[k].M.data), a.nops, a.basisH[k].size, time}")

for EL in ELlist:

    a.genHEBasis(k, subbasis, ET, EL)

    start = time.time()
    a.computeVhh(k, subbasis)
    end = time.time()

    # print("HE basis size", a.basisH[k].size)

    print("{",EL,",", len(a.Vhh[k].M.data),",",a.nops,",",a.basisH[k].size,",",
            end-start,"},")
