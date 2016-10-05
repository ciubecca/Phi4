import inspect
import os
import phi4
import sys
import math
import scipy
from statefuncs import *


def checkMatrix(matrix, basis, lookupbasis, fullmatrix, fullbasis, Emin, Emax, Vhh=False):

    Vred = matrix.todok()
    Vfull = fullmatrix.todok()

    helper = fullbasis.helper

    statePosRed = {}
    for i,state in enumerate(lookupbasis.stateList):
        statePosRed[tuple(helper.torepr2(state))] = i
        statePosRed[tuple(helper.torepr2(state)[::-1])] = i

    statePosFull = {}
    for I,state in enumerate(fullbasis.stateList):
        statePosFull[tuple(helper.torepr2(state))] = I
        statePosFull[tuple(helper.torepr2(state)[::-1])] = I


    stateListRed = [tuple(helper.torepr2(vi)) for vi in basis]
    stateListFull = [tuple(helper.torepr2(vJ)) for vJ in fullbasis]

    nnz = 0

    for i, vi in enumerate(basis):

        I = statePosFull[stateListRed[i]]

        for J,vJ in enumerate(fullbasis):

            # The matrix entry is not null. Therefore we can compare with the reduced
            # matrix
            if Emin<helper.energy(vJ)<=Emax and Vfull[I,J] != 0:

                try:
                    j = statePosRed[stateListFull[J]]
                except KeyError as err:
                    if Vhh==True:
                        continue
                    else:
                        raise err


                if abs(Vfull[I,J]-Vred[i,j]) > 10**(-10):
                    print("Numerical values don't match:")
                    print(Vfull[I,J], Vred[i,j])
                    print(vi, vJ)
                    print("Energies:", helper.energy(vi), helper.energy(vJ))
                    # raise ValueError

                # Count the number of non-zero entries
                nnz += 1

            else:
                try:
                    j = statePosRed[stateListFull[J]]
                    if Vred[i,j] != 0:
                        print("This entry should be zero")
                        print(vi, vJ)
                        print(helper.energy(vi), helper.energy(vJ))
                        raise ValueError
                except KeyError:
                    pass

    if nnz != matrix.nnz:
        print("Number of non-zero values don't match:")
        print(nnz, matrix.nnz)
        raise ValueError


k = 1
m = 1.

argv = sys.argv

args = "<L> <ET> <EL>"
if len(argv) < 4:
    print("python", argv[0], args)
    sys.exit(-1)

L = float(argv[1])
ET = float(argv[2])
EL = float(argv[3])

a = phi4.Phi4(m,L)

a.buildBasis(Emax=ET)

a.computePotential(k)

vset = [
[],
[(0, 2)],
# [(-1, 1), (1, 1)],
# [(-1, 1), (0, 2), (1, 1)],
# [(-2, 1), (-1, 1), (1, 1), (2, 1)]
]

# Build the reduced matrices VLh, Vhl
subbasis = Basis(k, vset, a.basis[k].helper)
print("subbasis size:", subbasis.size)

a.genHEBasis(k, subbasis, ET, EL)
print("HE basis size", a.basisH[k].size)

# Check that the total momentum of the states is 0
for v in a.basisH[k]:
    if sum(n*Zn for n,Zn in v) != 0:
        print("non-zero total momentum")
        raise ValueError

print("Emin, Emax = ", a.basisH[k].Emin, a.basisH[k].Emax)

eps = -1
print("Computing DeltaH")
a.computeDH2(k, subbasis, ET, EL, eps)

print("Computing Vhh")
a.computeVhh(k, subbasis)

# Build the full matrix up to cutoff ET
b = phi4.Phi4(m,L)
b.buildBasis(Emax=EL)
print("Full basis size:", b.basis[k].size)

b.computePotential(k)

tol = 10**-10

print("Checking Vhl")
checkMatrix(a.Vhl[k], subbasis, a.basisH[k], b.V[k][4].M, b.basis[k], ET, EL)

print("Checking VLh")
checkMatrix(a.VLh[k], a.basisH[k], a.basis[k], b.V[k][4].M, b.basis[k], 0-tol, ET-tol)


print("Checking Vhh")
checkMatrix(a.Vhh[k].M, a.basisH[k], a.basisH[k], b.V[k][4].M, b.basis[k],
        ET+tol, EL+tol, Vhh=True)
