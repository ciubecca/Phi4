import statefuncs
from  profile_support import *
import phi4
import renorm
import sys
import scipy
import math
import database

memdbg = False
# Whether the MonteCarlo integrals should be actually evaluated
test = False

if test:
    warnings.warn("Monte Carlo is OFF")
if memdbg:
    warnings.warn("Running with memory debugging")

# Whether we should save the results in the database data/spectra.db
saveondb = True
if not saveondb:
    warnings.warn("Saving on database is OFF")

m = 1
# Number of eigenvalues to compute per sector
neigs = 1
# List of parity quantum numbers
klist = (1,)
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

print("klist", klist)

print("maxntails", maxntails)

print("L, g, ET", L, g, ET)

EL = ratioELET*ET
ELp = ratioELpET*ET
ELpp = ratioELppELp*ELp

print("EL, ELp, ELpp", EL, ELp, ELpp)

# @profile
def main():

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

        # return

        # Select a set of tails and construct a Basis object
        vectorlist = [state for i,state in sorted(enumerate(a.basis[k]), key=lambda x:
                -abs(a.eigenvectors["raw"][k][0][x[0]]))]
        if maxntails != None:
            vectorlist = vectorlist[:maxntails]
        basisl = statefuncs.Basis(k, vectorlist, a.basis[k].helper)
        print("Total number of tails:", basisl.size)

        if memdbg:
            print("memory taken before computeLEVs", memory_usage())

        a.computeLEVs(k, basisl)


        if memdbg:
            print("memory taken before genHEBasis", memory_usage())

        print("Generating high energy basis...")
        # Generate the high-energy "selected" basis by passing a set of tails
        # and a maximum cutoff EL
        a.genHEBasis(k, EL=EL, ELp=ELp, ELpp=ELpp)
        print("Size of HE basis:", a.basisH[k].size)


        if memdbg:
            print("memory taken before computeHEVs", memory_usage())

        print("Computing high energy matrices...")
# Compute the matrices VLH, VHL, VHH, for the highest local cutoff ELmax.
# Later we will be varying EL, therefore taking submatrices of these.
# Computing VHH is expensive
        a.computeHEVs(k)

        a.computeEigval(k, ET, "renloc", eps=eps, neigs=neigs)
        print("Local ren vacuum:", a.eigenvalues["renloc"][k][0])
        eps = a.eigenvalues["renloc"][k][0]

        a.calcVV3([ELp], eps, test=test)

        a.computeEigval(k, ET, "rentails", EL=EL, ELp=ELp, ELpp=ELpp, eps=eps,
                neigs=neigs, memdbg=memdbg)
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
                    tailsComputedAtET=ET, maxntails=maxntails)
            db.insert(datadict=datadict, spec=a.eigenvalues["rentails"][k])



main()
