import statefuncs
import phi4
import renorm
import sys
import scipy
import math
import database

# Suppose we want to generate all the eigenvalues with L=10, g=1 with truncation energies
# ET = 10, 10.5, 11, ..., 20
# Then we should call this file as:
# eigsvsE.py 10 1 10 20

m = 1
# List of parity quantum numbers
klist = (-1,1)
minoverlap = 10**(-2)



argv = sys.argv
if len(argv) < 4:
    print(argv[0], "<L> <g> <ET>")
    sys.exit(-1)

L = float(argv[1])
g = float(argv[2])
ET = float(argv[3])


print("minoverlap:", minoverlap)


a = phi4.Phi4(m, L)
a.buildBasis(Emax=ET)


for k in klist:

    # Compute the potential matrices in the low-energy space below ET
    a.computePotential(k)

    print("k=", k)
    print("Full basis size: ", a.basis[k].size)

    a.setCouplings(g4=g)
    print("g=", g)


    print("Computing raw eigenvalues for highest cutoff")
    a.computeEigval(k, ET, "raw", neigs=10)


    # Select a set of tails and construct a Basis object
    vectorlist = [state for i,state in enumerate(a.basis[k])
        if abs(a.eigenvectors["raw"][k][0][i]) > minoverlap]
    basisl = statefuncs.Basis(k, vectorlist, a.basis[k].helper)
    print("Total number of tails:", basisl.size)

    a.genHEBasis(k,basisl,ET+1)

    print("Computing low energy matrices...")
# Compute the matrices VLH, VHL, VHH, for the highest local cutoff ELmax.
# Later we will be varying EL, therefore taking submatrices of these.
# Computing VHH is expensive
    a.computeLEVs(k)

    print(a.V2V4[k].M.todense())
