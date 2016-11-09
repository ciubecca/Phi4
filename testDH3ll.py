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


def checkMatrix(matrix, basis, lookupbasis, fullmatrix, fullbasis, Emin, Emax, Vll=False):

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
                    if Vll==True:
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

args = "<L> <ET> <EL> <ELp> <ELpmax>"
if len(argv) < 6:
    print("python", argv[0], args)
    sys.exit(-1)

L = float(argv[1])
ET = float(argv[2])
EL = float(argv[3])
ELp = float(argv[4])
ELpmax = float(argv[5])

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

a.genHEBases(k, basisl, EL=EL, ELp=ELpmax)
print("HE basis size", a.basisH[k].size)

a.computeLEVs(k)

# Check that the total momentum of the states is 0
for v in a.basisH[k]:
    if sum(n*Zn for n,Zn in v) != 0:
        print("non-zero total momentum")
        raise ValueError

print("Emin, Emax = ", a.basisH[k].Emin, a.basisH[k].Emax)

print("Computing HE matrices")
a.computeHEVs(k)


# Build the full matrix up to cutoff ET
b = phi4.Phi4(m,L)
b.buildBasis(Emax=EL)
print("Full basis size:", b.basis[k].size)

b.computePotential(k)

tol = 10**-10


# Subset of the selected low energy states
basisl = a.basisl[k]
helper = basisl.helper
vectorlist = [v for v in basisl if helper.energy(v)<=ET]
subbasisl = Basis(k, vectorlist, helper)
helper = a.basisH[k].helper
subbasish = a.basish[k].sub(lambda v: ET<helper.energy(v)<=ELp)

Vhl = {}
for n in (2,4):
    Vhl[n] = a.Vhl[n][k].sub(subbasisl,subbasish).M

for n in (2,4):
    print("Checking Vhl {}".format(n))
    checkMatrix(Vhl[n], subbasisl, subbasish, b.V[k][n].M, b.basis[k], ET, ELp)



# Subsets of the selected high energy states
helper = a.basisH[k].helper
subbasisH = a.basisH[k].sub(lambda v: ET<helper.energy(v)<=EL)

VHl = a.VHl[k].sub(subbasisl, subbasisH).M

print("Checking VHl")
checkMatrix(VHl, subbasisl, subbasisH, b.V[k][4].M, b.basis[k], ET, EL)




helper = a.basis[k].helper
# Subset of the full low energy states
subbasisL = a.basis[k].sub(lambda v: helper.energy(v)<=ET)

VLH = a.VLH[k].sub(subbasisH, subbasisL).M

print("Checking VLH")
checkMatrix(VLH, subbasisH, a.basis[k], b.V[k][4].M, b.basis[k], 0-tol, ET)



# print("Checking Vhh")
# checkMatrix(a.Vhh[k].M, a.basisH[k], a.basisH[k], b.V[k][4].M, b.basis[k], EL+tol, Vhh=True)


print("Building DH3ll")

loc3mix = False
nonloc3mix = False
loc3 = False
eps = -1

DH3ll = a.computeDH3(subbasisl, k, ET, ELp, ELpp=None, eps=eps, loc3=False, loc3mix=False,
        nonloc3mix=False)

basis = b.basis[k]
energyArr = array(basis.energyList)
propagator = scipy.sparse.spdiags(1/(eps-energyArr), 0, basis.size, basis.size)
proj = scipy.sparse.spdiags(array([int(ET<e<ELp) for e in energyArr]), 0, basis.size,
        basis.size)
Vfull = b.V[k][4].M
DH3Full = Vfull*propagator*proj*Vfull*propagator*proj*Vfull*g**3


print("Checking DH3ll")
checkMatrix(DH3ll, subbasisl, subbasisl, DH3Full, basis, 0-tol, ET, Vll=True)
