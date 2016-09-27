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


# Build the full matrix up to cutoff ET
b = phi4.Phi4(m,L)
b.buildBasis(Emax=Ebar, occmax=occmax)
b.buildMatrix()

# Index list of basis elements in the selected low energy basis
subIndexList = [i for i in range(b.basis[k].size) if b.basis[k].stateList[i] in vset]

# Full matrix up to cutoff ET
Vfull = b.V[k][4].M.todok()

# Reduced low-high energy matrix
Vred = a.VLH[k].M.todok()


# Get all the indices j such that V[i,j]!=0 in the full matrix
HEindexSet = set()
for i in subIndexList:
    for j,v in enumerate(b.basis[k].stateList):
        if b.basis[k].energyList[j] > Emax and Vfull[i,j] != 0:
            HEindexSet.add(j)


# Get all the indices j in the full basis which correspond to states
# in the reduced high energy basis
HEindexSet2 = set()
for i, state in enumerate(a.basisH[k]):
    try:
        HEindexSet2.add(b.basis[k].lookup(state))
    except KeyError:
        print("KeyError:", state)
        # print(a.basisH[k].energyList[i])



# HEindexList2 = []
# d = {}
# for j, state in enumerate(a.basisH[k]):
    # i = b.basis[k].lookup(state)
    # if i in HEindexList2:
        # print("Already present:", d[i], j)
        # print("Already present:", a.basisH[k].stateList[d[i]], b.basis[k].stateList[i], state)
    # d[i] = j
    # HEindexList2.append(i)


# Check that the number of reduced high energy states corresponds to the number
# of non-zero matrix elements of the full matrix on the reduced low energy subspace
print(len(HEindexSet)-len(HEindexSet2))


# Check that the numerical values of the reduced and full matrices correspond
for i, vl in enumerate(vset):
    for j, vh in enumerate(a.basisH[k]):
        i2 = b.basis[k].lookup(vl)
        j2 = b.basis[k].lookup(vh)

        if abs(Vred[i,j]-Vfull[i2,j2]) > 10**(-10):
            print("Difference:", i, j)
            print("Values:", Vred[i,j], Vfull[i2,j2])
            print("States:", vl, vh)
            print("Energies:", b.basis[k].energyList[i2],b.basis[k].energyList[j2])

