import statefuncs
import phi4
import renorm
import sys
import scipy
import math
import database


# Whether we should save the results in the database data/spectra.db
saveondb = True
# saveondb = False
m = 1
# Number of eigenvalues to compute per sector
neigs = 1
# List of parity quantum numbers
klist = (-1,1)
# Minimum overlap with the raw vacuum for selecting a state in the tails
minoverlap = 10**(-4)


# How EL depends on ET
def el(ET):
    return ET*2.5

argv = sys.argv
if len(argv) < 4:
    print(argv[0], "<L> <g> <ET>")
    sys.exit(-1)

L = float(argv[1])
g = float(argv[2])
ET = float(argv[3])


if saveondb:
    db = database.Database()

a = phi4.Phi4(m, L)
a.buildBasis(Emax=ET)


for k in klist:

    # Compute the potential matrices in the low-energy space below ET
    a.computePotential(k)

    print("k=", k)
    print("Full basis size: ", a.basis[k].size)

    a.setCouplings(g4=g)
    print("g=", g)


# Compute the raw eigenvalues for cutoff ET
    a.computeEigval(k, ET, "raw", neigs=neigs)
    print("Raw vacuum:", a.eigenvalues["raw"][k][0])
    eps = a.eigenvalues["raw"][k][0]



    # Select a set of tails and construct a Basis object
    vectorlist = [state for i,state in enumerate(a.basis[k])
        if abs(a.eigenvectors["raw"][k][0][i]) > minoverlap]
    basisl = statefuncs.Basis(k, vectorlist, a.basis[k].helper)
    print("Total number of tails:", basisl.size)


    EL = el(ET)

    print("ET={}, EL={}".format(ET,EL))

    print("Generating high energy basis...")
    # Generate the high-energy "selected" basis by passing a set of tails
    # and a maximum cutoff EL
    a.genHEBasis(k, basisl, EL)
    print("Size of HE basis:", a.basisH[k].size)

    print("Computing high energy matrices...")
# Compute the matrices VLH, VHL, VHH, for the highest local cutoff ELmax.
# Later we will be varying EL, therefore taking submatrices of these.
# Computing VHH is expensive
    a.computeHEVs(k, EL)


    # Checks if the eigenvalues we are going to compute are alredy present in
#the database. If yes, skip them
    if saveondb:
        approxQuery = {"g":g, "L":L, "ET":ET, "EL":EL}
        exactQuery = {"k":k, "ren":"rentails"}
        if db.getObjList('spec', approxQuery=approxQuery,
                exactQuery=exactQuery) != []:
            print("Eigenvalues already present")
            continue


    if saveondb:
        db.insert(k=k, ET=ET, L=L, ren="raw", g=g,
                spec=a.eigenvalues["raw"][k], eigv=a.eigenvectors["raw"][k],
                basisSize=a.compSize, neigs=neigs)


# Compute "local" renormalized eigenvalues for cutoff ET
# Since we are passing EL=ET to the method call, the matrices VHL, VHH will be computed
# only in the local approximation
    a.computeEigval(k, ET, "renloc", ET, eps, neigs=neigs)
    print("Local ren vacuum:", a.eigenvalues["renloc"][k][0])
    eps = a.eigenvalues["renloc"][k][0]

    if saveondb:
        db.insert(k=k, ET=ET, L=L, ren="renloc", g=g,
                spec=a.eigenvalues["renloc"][k], eigv=a.eigenvectors["renloc"][k],
                basisSize=a.compSize, neigs=neigs, EL=ET, eps=eps)

# Compute renormalized eigenvalues by computing the fully "non-local" corrections
# to VHL, VHH up to cutoff EL
    a.computeEigval(k, ET, "rentails", EL, eps, neigs=neigs)
    print("Non-Local ren vacuum:", a.eigenvalues["rentails"][k][0])

    print("Number of tails:", a.ntails)

    if saveondb:
        db.insert(k=k, ET=ET, L=L, ren="rentails", g=g,
                spec=a.eigenvalues["rentails"][k], eigv=a.eigenvectors["rentails"][k],
                basisSize=a.compSize, neigs=neigs, EL=EL, ntails=a.ntails, eps=eps)
