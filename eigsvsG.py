import statefuncs
from  profile_support import *
import phi4
import renorm
import sys
import scipy
import math
import database

glist = scipy.linspace(0.5, 3, 6)
print("glist", glist)

memdbg = True
if memdbg:
    warnings.warn("Running with memory debugging")
# Whether the MonteCarlo integrals should be actually evaluated
test = False
if test:
    warnings.warn("Monte Carlo is OFF")
loc3 = False
if not loc3:
    warnings.warn("Not including local correction to DH3")


# Whether we should save the results in the database data/spectra.db
saveondb = True
if not saveondb:
    warnings.warn("Saving on database is OFF")

m = 1
# Number of eigenvalues to compute per sector
neigs = 1
# List of parity quantum numbers
klist = (1,)

# Ratio between EL and ET
ratioELET = 3
# Ratio between ELp and ET
ratioELpET = 2
# Ratio between ELpp and ELp
ratioELppELp = 1.5


argv = sys.argv
if len(argv) < 3:
    print(argv[0], "<L> <ET>")
    sys.exit(-1)

L = float(argv[1])
ET = float(argv[2])

print("klist", klist)

print("L, ET", L, ET)

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

        a.setglist(glist=glist)


# Compute the raw eigenvalues for cutoff ET
        a.computeEigval(k, ET, "raw", neigs=neigs)
        epsraw = {g: a.eigenvalues[g]["raw"][k][0] for g in glist}
        print("Raw vacuum:", epsraw)

        # return

        # Always consider maximal set of tails
        basisl = a.basis[k]
        print("Total number of tails:", basisl.size)

        if memdbg:
            print("memory taken before computeLEVs", memory_usage())

        a.computeLEVs(k, basisl, loc3=loc3)

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

        a.computeEigval(k, ET, "renloc", eps=epsraw, neigs=neigs)
        eps = {g: a.eigenvalues[g]["renloc"][k][0] for g in glist}
        print("Local ren vacuum:", eps)

        if loc3:
            a.calcVV3(ELp, eps, test=test)

        a.computeEigval(k, ET, "rentails", EL=EL, ELp=ELp, ELpp=ELpp, eps=eps,
                neigs=neigs, memdbg=memdbg, loc3=loc3)
        print("Non-Local ren vacuum:", {g: a.eigenvalues[g]["rentails"][k][0]
            for g in glist})

        print("Number of tails:", a.ntails)

        if saveondb:
            for g in glist:
                datadict = dict(k=k, ET=ET, L=L, ren="raw", g=g, neigs=neigs,
                        basisSize=a.compSize)
                db.insert(datadict=datadict, spec=a.eigenvalues[g]["raw"][k])

                datadict = dict(k=k, ET=ET, L=L, ren="renloc", g=g, eps=epsraw[g],
                        neigs=neigs, basisSize=a.compSize)
                db.insert(datadict=datadict, spec=a.eigenvalues[g]["renloc"][k])

                datadict = dict(k=k, ET=ET, L=L, ren="rentails", g=g, EL=EL, ELp=ELp,
                        ELpp=ELpp, ntails=a.ntails, eps=eps[g], neigs=neigs,
                        basisSize=a.compSize, tailsComputedAtET=ET, maxntails=None)
                db.insert(datadict=datadict, spec=a.eigenvalues[g]["rentails"][k])



main()
