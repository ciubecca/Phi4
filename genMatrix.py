import phi4
import sys
from sys import getsizeof as size
import math
import scipy
from statefuncs import *

k = 1
m = 1.

ntails = 200

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

print("basis size :", a.basis[k].size)

a.computePotential(k)

basisl = Basis(k, a.basis[k].stateList[:ntails], a.basis[k].helper)

print("Number of tails:", basisl.size)

a.genHEBasis(k, basisl, EL)

print("HE basis size", a.basisH[k].size)

print("Computing high energy matrices...")
a.computeHEVs(k, EL)

# VhhList = a.VhhHalfList[k]
# print("Size of Vhh in memory:",
    # sum(int((size(Vhh.data)+size(Vhh.indptr)+size(Vhh.indices))) for Vhh in VhhList))

print("Computing DeltaH...")
eps = -1
a.setCouplings(g4=1)
a.computeDeltaH(k, ET, EL, eps, 200)

# a.saveMatrix(k=k, Emax=Emax)
