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

nmax =  a.helper.nmax

a.buildMatrix(k)


vset = [
[],
[(0, 2)],
[(-1, 1), (1, 1)],
[(-1, 1), (0, 2), (1, 1)],
[(-2, 1), (-1, 1), (1, 1), (2, 1)]
]

subbasis = Basis(k, vset, a.helper)
print("subbasis size:", subbasis.size)

a.computeDH2(k, subbasis, Emax, Ebar)

print("HE basis size", a.basisH[k].size)


# Build the full matrix up to cutoff ET
b = phi4.Phi4(m,L)
b.buildBasis(Emax=Ebar, occmax=occmax)
b.buildMatrix(k)

# Full matrix up to cutoff ET
Vfull = b.V[k][4].M.todok()

# Reduced low-high energy matrix
Vred = a.VLH[k].M.todok()

# Check that the numerical values of the reduced and full matrices correspond
for i, vl in enumerate(vset):
    for j, vh in enumerate(a.basisH[k]):
        i2 = b.basis[k].lookup(vl)
        j2 = b.basis[k].lookup(vh)

        if abs(Vred[i,j]-Vfull[i2,j2]) > 10**(-10):
            print("VLH difference:", i, j)
            print("Values:", Vred[i,j], Vfull[i2,j2])
            print("States:", vl, vh)
            print("Energies:", b.basis[k].energyList[i2],b.basis[k].energyList[j2])
            raise ValueError


# Reduced high-low energy matrix
Vred = a.VHL[k].M.todok()

# Check that the numerical values of the reduced and full matrices correspond
for i, vh in enumerate(a.basisH[k]):
    for j, vl in enumerate(a.basis[k]):
        i2 = b.basis[k].lookup(vh)
        j2 = b.basis[k].lookup(vl)

        if abs(Vred[i,j]-Vfull[i2,j2]) > 10**(-10):
            print("VHL difference:", i, j)
            print("Values:", Vred[i,j], Vfull[i2,j2])
            print("States:", vh, vl)
            print("Energies:", b.basis[k].energyList[i2],b.basis[k].energyList[j2])
            raise ValueError

