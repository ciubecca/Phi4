import inspect
import os
import phi4
import sys
import math
import scipy
from statefuncs import *


k = +1
m = 1.
g = .1


minoverlap = 10**-3


def checkMatrix(matrix, basis, lookupbasis, fullmatrix, fullbasis, Emax, Vhh=False):

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
            if helper.energy(vJ)<=Emax and Vfull[I,J] != 0:

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
                    raise ValueError

                # Count the number of non-zero entries
                nnz += 1

            else:
                try:
                    j = statePosRed[stateListFull[J]]
                    if Vred[i,j] != 0:
                        print("This entry should be zero")
                        print(Vfull[I,J], Vred[i,j])
                        print(vi, vJ)
                        print(helper.energy(vi), helper.energy(vJ))
                        raise ValueError
                except KeyError:
                    pass

    if nnz != matrix.nnz:
        print("Number of non-zero values don't match:")
        print(nnz, matrix.nnz)
        raise ValueError



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

a.setCouplings(0,0,g)
a.computeEigval(k, ET, "raw")

vectorlist = [state for i,state in enumerate(a.basis[k])
        if abs(a.eigenvectors["raw"][1][0][i]) > minoverlap]
print(sorted(occn(state) for state in vectorlist))
basisl = Basis(k, vectorlist, a.basis[k].helper)
print("subbasis size:", basisl.size)

a.genHEBasis(k, basisl, EL)
print("HE basis size", a.basisH[k].size)

# Check that the total momentum of the states is 0
for v in a.basisH[k]:
    if sum(n*Zn for n,Zn in v) != 0:
        print("non-zero total momentum")
        raise ValueError

print("Emin, Emax = ", a.basisH[k].Emin, a.basisH[k].Emax)

print("Computing HE matrices")
a.computeHEVs(k, EL)


# Build the full matrix up to cutoff ET
b = phi4.Phi4(m,L)
b.buildBasis(Emax=EL)
print("Full basis size:", b.basis[k].size)

b.computePotential(k)

tol = 10**-10

print("Checking Vhl")
checkMatrix(a.Vhl[k].M, basisl, a.basisH[k], b.V[k][4].M, b.basis[k], EL)

print("Checking VLh")
checkMatrix(a.VLh[k].M, a.basisH[k], a.basis[k], b.V[k][4].M, b.basis[k], ET+tol)


print("Checking Vhh")
checkMatrix(a.Vhh[k].M, a.basisH[k], a.basisH[k], b.V[k][4].M, b.basis[k], EL+tol, Vhh=True)
