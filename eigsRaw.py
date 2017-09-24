import statefuncs
from  profile_support import *
import phi4
import renorm
import sys
import scipy
import numpy as np
import math
import database

memdbg = False
if memdbg:
    warnings.warn("Running with memory debugging")
# Whether the MonteCarlo integrals should be actually evaluated
test = False
if test:
    warnings.warn("Monte Carlo is OFF")
loc3 = True
if not loc3:
    warnings.warn("Not including local correction to DH3")


# Whether we should save the results in the database data/spectra.db
saveondb = False
if not saveondb:
    warnings.warn("Saving on database is OFF")

m = 1
# Number of eigenvalues to compute per sector
neigs = 6

# Ratio between EL and ET
ratioELET = 3
# Ratio between ELp and ET
ratioELpET = 2
# Ratio between ELpp and ELp
ratioELppELp = 1.5


argv = sys.argv
if len(argv) < 5:
    print(argv[0], "<k> <L> <ET> <theta>")
    sys.exit(-1)

k = int(argv[1])
L = float(argv[2])
ET = float(argv[3])
theta = float(argv[4])


glist = [0.5*np.exp(1J*theta)]

print(glist)


print("k, L, ET", k, L, ET)

EL = ratioELET*ET
ELp = ratioELpET*ET
ELpp = ratioELppELp*ELp

print("EL, ELp, ELpp", EL, ELp, ELpp)

def main():

    if saveondb:
        db = database.Database()

    a = phi4.Phi4(m, L, k)
    a.buildBasis(Emax=ET)

    # Compute the potential matrices in the low-energy space below ET
    a.computePotential()

    print("Full basis size: ", a.basis.size)

    a.setglist(glist=glist)

# Compute the raw eigenvalues for cutoff ET
    a.computeEigval(ET, "raw", neigs=neigs)
    epsraw = {g: a.eigenvalues[g]["raw"][0] for g in glist}
    print("Raw vacuum:", epsraw)


    if saveondb:
        for g in glist:
            datadict = dict(k=k, ET=ET, L=L, ren="raw", g=g, neigs=neigs,
                    basisSize=a.compSize, finiteL=True)
            db.insert(datadict=datadict, spec=a.eigenvalues[g]["raw"])

main()
