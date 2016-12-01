import statefuncs
import phi4
import renorm
import sys
import scipy
import math
import gc
import database
from statefuncs import *


m = 1
# List of parity quantum numbers
k = 1

# Maximum number of tails (in order of overlap with the vacuum) to include
maxntails = 100
ntails = 50

# Ratio between ELpp and ELp
ratio3 = 1.1
# Ratio between EL and ET
# ratio2 = 2
neigs = 10

loc3 = True

argv = sys.argv
if len(argv) < 6:
    print(argv[0], "<L> <g> <ET> <EL> <ELp>")
    sys.exit(-1)

L = float(argv[1])
g = float(argv[2])
ET = float(argv[3])
EL = float(argv[4])
ELp = float(argv[5])

ELpp = ratio3*ELp

print("maxntails:", maxntails)
print("ELpp/ELp:", ratio3)


compDH3 = {}
compDH2 = {}
vectorList = {}
Vllcomp = {}
V0V4comp = {}
V2V4comp = {}
V4V4comp = {}

for compntails in (ntails,maxntails):
    print("computing DH3 at most for ntails=", compntails)

    a = phi4.Phi4(m, L)
    a.buildBasis(Emax=ET)

    Vllcomp[compntails] = {}


    print("Full basis size:", a.basis[k].size)

    # Compute the potential matrices in the low-energy space below ET
    a.computePotential(k)

    a.setCouplings(g4=g)

    a.computeEigval(k, ET, "raw")
    eps = a.eigenvalues["raw"][k][0]

    # Select a set of tails and construct a Basis object, ordered in overlap with
    # the vacuum
    vectorlist = [state for i,state in sorted(enumerate(a.basis[k]), key=lambda x:
            -abs(a.eigenvectors["raw"][k][0][x[0]]))][:compntails]
    basisl = statefuncs.Basis(k, vectorlist, a.basis[k].helper, orderEnergy=False)
    # print(vectorlist)
    vectorList[compntails] = basisl.stateList
    print("basisl size:", basisl.size)

    # Generate the high-energy "selected" basis by passing a set of tails
    # and a maximum cutoff EL
    a.genHEBases(k, basisl, EL=EL, ELpp=ELpp)
    print("HE basis sizes:", a.basisH[k].size, a.basish[k].size)

    a.computeLEVs(k)

# Compute the matrices VLH, VHL, VHH, for the highest local cutoff ELmax.
# Later we will be varying EL, therefore taking submatrices of these.
# Computing VHH is expensive
    a.computeHEVs(k)


# Compute "local" renormalized eigenvalues for cutoff ET
# Since we are passing EL=ET to the method call, the matrices VHL, VHH will be computed
# only in the local approximation
    a.computeEigval(k, ET, "renloc", EL=ET, eps=eps)
    eps = a.eigenvalues["renloc"][k][0]

    if loc3:
        a.calcVV3([ELp], eps)

    if compntails==maxntails:
        a.precomputeDH3(k, ET, ELp, ELpp, eps, loc3=loc3)

    a.computeEigval(k, ET, "rentails", EL=EL, ELp=ELp, ELpp=ELpp, eps=eps, neigs=neigs,
                maxntails=ntails, loc3=loc3)
    print("Non-Local ren vacuum:", a.eigenvalues["rentails"][k][0])
    print("Number of tails:", a.ntails)

    compDH3[compntails] = a.compDH3[k]
    compDH2[compntails] = a.compDH2[k]
    for n in (0,2,4,6):
        Vllcomp[compntails][n] = a.Vllcomp[k][n]


    V0V4comp[compntails] = a.V0V4comp[k]
    V2V4comp[compntails] = a.V2V4comp[k]
    V4V4comp[compntails] = a.V4V4comp[k]



print("DH3")
print(compDH3[ntails][k]-compDH3[maxntails][k])

# print("V0V4")
# print(V0V4comp[ntails]-V0V4comp[maxntails])
# print("V2V4")
# print(V2V4comp[ntails]-V2V4comp[maxntails])
# print("V4V4")
# print(V4V4comp[ntails]-V4V4comp[maxntails])
# for n in (0,2,4,6):
    # print("Vll", n)
    # print(Vllcomp[ntails][n]-Vllcomp[maxntails][n])

