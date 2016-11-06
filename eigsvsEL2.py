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

# Whether we should save the results in the database data/spectra.db
saveondb = True
# saveondb = False
m = 1
# Number of eigenvalues to compute per sector
neigs = 1
# List of parity quantum numbers
klist = (-1,1)
# Minimum overlap with the raw vacuum for selecting a state in the tails
minoverlap = 10**(-2)



argv = sys.argv
if len(argv) < 6:
    print(argv[0], "<L> <g> <ET> <EL2min> <EL2max>")
    sys.exit(-1)

L = float(argv[1])
g = float(argv[2])
ET = float(argv[3])
EL2min = float(argv[4])
EL2max = float(argv[5])
# ELETdiff = float(argv[5])

EL3 = EL2min

EL2list = scipy.linspace(EL2min, EL2max, (EL2max-EL2min)*2+1)
print("EL2list:", EL2list)

print("minoverlap:", minoverlap)

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


    print("Computing raw eigenvalues for highest cutoff")
    a.computeEigval(k, ET, "raw", neigs=neigs)


    # Select a set of tails and construct a Basis object
    vectorlist = [state for i,state in enumerate(a.basis[k])
        if abs(a.eigenvectors["raw"][k][0][i]) > minoverlap]
    basisl = statefuncs.Basis(k, vectorlist, a.basis[k].helper)
    print("Total number of tails:", basisl.size)


    print("Generating high energy basis...")
    # Generate the high-energy "selected" basis by passing a set of tails
    # and a maximum cutoff EL
    a.genHEBasis(k, basisl, EL2max)
    print("Size of HE basis:", a.basisH[k].size)


    a.computeLEVs(k)


    print("Computing high energy matrices...")
# Compute the matrices VLH, VHL, VHH, for the highest local cutoff ELmax.
# Later we will be varying EL, therefore taking submatrices of these.
# Computing VHH is expensive
    a.computeHEVs(k, EL2max)


    # Checks if the eigenvalues we are going to compute are alredy present in
#the database. If yes, skip them
    # if saveondb:
        # approxQuery = {"g":g, "L":L, "ET":ET, "EL":EL}
        # exactQuery = {"k":k, "ren":"rentails"}
        # if db.getObjList('spec', approxQuery=approxQuery,
                # exactQuery=exactQuery) != []:
            # print("Eigenvalues already present")
            # continue

# Compute the raw eigenvalues for cutoff ET
    a.computeEigval(k, ET, "raw", neigs=neigs)
    print("Raw vacuum:", a.eigenvalues["raw"][k][0])
    eps = a.eigenvalues["raw"][k][0]

    if saveondb:
        db.insert(k=k, ET=ET, L=L, ren="raw", g=g,
                spec=a.eigenvalues["raw"][k],
                basisSize=a.compSize, neigs=neigs)


# Compute "local" renormalized eigenvalues for cutoff ET
# Since we are passing EL=ET to the method call, the matrices VHL, VHH will be computed
# only in the local approximation
    a.computeEigval(k, ET, "renloc", EL=ET, eps=eps, neigs=neigs)
    print("Local ren vacuum:", a.eigenvalues["renloc"][k][0])
    eps = a.eigenvalues["renloc"][k][0]

    if saveondb:
        db.insert(k=k, ET=ET, L=L, ren="renloc", g=g,
                spec=a.eigenvalues["renloc"][k],
                basisSize=a.compSize, neigs=neigs, EL=ET, eps=eps)


    # a.calcVV3(ETlist=[EL3], eps=eps)


    for EL in EL2list:

        # EL = ELmax
        print("EL={}".format(EL))


# Compute renormalized eigenvalues by computing the fully "non-local" corrections
# to VHL, VHH up to cutoff EL
        a.computeEigval(k, ET, "rentails", EL, eps=eps, neigs=neigs, EL3=EL3)
        print("Non-Local ren vacuum:", a.eigenvalues["rentails"][k][0])

        print("Number of tails:", a.ntails)

        if saveondb:
            db.insert(k=k, ET=ET, L=L, ren="rentails", g=g, minoverlap=minoverlap,
                    spec=a.eigenvalues["rentails"][k], EL3=EL3,
                    basisSize=a.compSize, neigs=neigs, EL=EL, ntails=a.ntails, eps=eps)


    del a.VLh[k]
