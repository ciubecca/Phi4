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
print("nmax", nmax)


# print(a.basis[k])
print("basis size :", a.basis[k].size)

# print(a.basis[1][-1])
# print(a.basis[1].stateDlists[2][-1])

a.buildMatrix()
# print(a.basis[k])

subbasis = a.basis[k].sub(lambda v: occn(v)<=4)
print("subbasis size:", subbasis.size)

a.computeDH2(k, subbasis, Emax, Ebar)


# print(a.basisH[1].energyList)
print("HE basis size", a.basisH[k].size)
# print("Emax", a.basisH[1].Emin)
# print("Emin", a.basisH[1].Emax)

print(a.basisH[k].stateList[-10:])
print(a.basisH[k].energyList[-10:])

# print("Estimated memory use of basis:", 16*(2*nmax+1)*a.basisH[1].size/(10^6), "M")

# a.saveMatrix(k=k)

# print(a.V[k][4].M.todense())

# print(a.V[k][4].M.todense())
