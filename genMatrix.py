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

args = "<L> <Emax> <Ebar> <?occmax>"
if len(argv) < 4:
    print("python", argv[0], args)
    sys.exit(-1)

L = float(argv[1])
Emax = float(argv[2])
Ebar = float(argv[3])
try:
    occmax = int(argv[4])
except IndexError:
    occmax = None

a = phi4.Phi4(m,L)

a.buildBasis(Emax=Emax, occmax=occmax)

# print(a.basis[k])
print("basis size :", a.basis[k].size)

a.buildMatrix()

vset = [
[],
[(0, 2)],
[(-1, 1), (1, 1)],
[(-1, 1), (0, 2), (1, 1)],
[(-2, 1), (-1, 1), (1, 1), (2, 1)]
]


subbasis = Basis(k, vset, a.helper, gendlist=True)
print("subbasis size:", subbasis.size)

a.computeDH2(k, subbasis, Emax, Ebar)


print("HE basis size", a.basisH[k].size)


a.saveMatrix(k=k, Emax=Emax)
