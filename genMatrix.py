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

args = "<L> <ET> <EL> <?occmax>"
if len(argv) < 4:
    print("python", argv[0], args)
    sys.exit(-1)

L = float(argv[1])
ET = float(argv[2])
EL = float(argv[3])
try:
    occmax = int(argv[4])
except IndexError:
    occmax = None

a = phi4.Phi4(m,L)

a.buildBasis(Emax=ET, occmax=occmax)

print("basis size :", a.basis[k].size)

a.computePotential(k)

# vset = [
# [],
# [(0, 2)],
# [(-1, 1), (1, 1)],
# [(-1, 1), (0, 2), (1, 1)],
# [(-2, 1), (-1, 1), (1, 1), (2, 1)]
# ]
# subbasis = Basis(k, vset, a.helper, gendlist=True)


subbasis = Basis(k, a.basis[k].stateList[:50], a.basis[k].helper)

a.genHEBasis(k, subbasis, ET, EL)

print("HE basis size", a.basisH[k].size)

eps = -1
a.computeDH2(k, subbasis, ET, EL, eps)

# a.computeVHH(k, subbasis, Emax, Ebar)

# a.saveMatrix(k=k, Emax=Emax)
