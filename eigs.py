import statefuncs
import phi4
import renorm
import sys
import scipy
import math
import database

loc3 = False

# Whether we should save the results in the database data/spectra.db
saveondb = False
# saveondb = False
m = 1
# Number of eigenvalues to compute per sector
neigs = 1
# List of parity quantum numbers
klist = (1,)
# Minimum overlap with the raw vacuum for selecting a state in the tails
maxntails = None

# Ratio between EL and ET
ratioELET = 3
# Ratio between ELp and ET
ratioELpET = 2
# Ratio between ELpp and ELp
ratioELppELp = 1.5

argv = sys.argv
if len(argv) < 4:
    print(argv[0], "<L> <g> <ET>")
    sys.exit(-1)

L = float(argv[1])
g = float(argv[2])
ET = float(argv[3])

EL = ratioELET*ET
ELp = ratioELpET*ET
ELpp = ratioELppELp*ELp

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


    # Select a set of tails and construct a Basis object, ordered in overlap with
    # the vacuum
    vectorlist = [state for i,state in sorted(enumerate(a.basis[k]), key=lambda x:
            -abs(a.eigenvectors["raw"][k][0][x[0]]))]
    if maxntails!=None:
        vectorlist = vectorlist[:maxntails]
    basisl = statefuncs.Basis(k, vectorlist, a.basis[k].helper, orderEnergy=False)
    print("ntails:", basisl.size)


    print("Generating high energy basis...")
    # Generate the high-energy "selected" basis by passing a set of tails
    # and a maximum cutoff EL
    a.genHEBases(k, basisl, EL=EL, ELpp=ELpp)
    print("Size of HE basis:", a.basisH[k].size)

    print("Computing high energy matrices...")
    a.computeHEVs(k)
# Compute the matrices VLH, VHL, VHH, for the highest local cutoff ELmax.
# Later we will be varying EL, therefore taking submatrices of these.
# Computing VHH is expensive


    a.computeLEVs(k)


# Compute "local" renormalized eigenvalues for cutoff ET
# Since we are passing EL=ET to the method call, the matrices VHL, VHH will be computed
# only in the local approximation
    a.computeEigval(k, ET, "renloc", eps=eps, neigs=neigs)
    print("Local ren vacuum:", a.eigenvalues["renloc"][k][0])
    eps = a.eigenvalues["renloc"][k][0]

    if loc3:
        a.calcVV3([ELp], eps)

# Compute renormalized eigenvalues by computing the fully "non-local" corrections
# to VHL, VHH up to cutoff EL
    a.computeEigval(k, ET, "rentails", EL=EL, ELp=ELp, ELpp=ELpp, eps=eps, neigs=neigs,
            loc3=loc3)
    print("Non-Local ren vacuum:", a.eigenvalues["rentails"][k][0])

    print("Number of tails:", a.ntails)

    if saveondb:
        datadict = dict(k=k, ET=ET, L=L, ren="raw", g=g, neigs=neigs,
                basisSize=a.compSize)
        db.insert(datadict=datadict, spec=a.eigenvalues["raw"][k])


        datadict = dict(k=k, ET=ET, L=L, ren="renloc", g=g, eps=eps, neigs=neigs,
                basisSize=a.compSize)
        db.insert(datadict=datadict, spec=a.eigenvalues["renloc"][k])


        datadict = dict(k=k, ET=ET, L=L, ren="rentails", g=g, EL=EL, ELp=ELp, ELpp=ELpp,
                ntails=a.ntails, eps=eps, neigs=neigs, basisSize=a.compSize,
                tailsComputedAtET=ETmax, maxntails=maxntails)
        db.insert(datadict=datadict, spec=a.eigenvalues["rentails"][k])


