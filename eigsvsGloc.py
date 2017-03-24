import statefuncs
from  profile_support import *
import phi4
import renorm
import sys
import scipy
import math
import database

glist = scipy.linspace(1, 3, 3)
print("glist", glist)

memdbg = False
if memdbg:
    warnings.warn("Running with memory debugging")

# Whether we should save the results in the database data/spectra.db
saveondb = True
if not saveondb:
    warnings.warn("Saving on database is OFF")

m = 1
# Number of eigenvalues to compute per sector
neigs = 6

argv = sys.argv
if len(argv) < 4:
    print(argv[0], "<k> <L> <ET>")
    sys.exit(-1)

k = int(argv[1])
L = float(argv[2])
ET = float(argv[3])


print("k, L, ET", k, L, ET)


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


    a.computeEigval(ET, "renloc", eps=epsraw, neigs=neigs)
    eps = {g: a.eigenvalues[g]["renloc"][0] for g in glist}
    print("Local ren vacuum:", eps)

    if saveondb:
        for g in glist:
            datadict = dict(k=k, ET=ET, L=L, ren="raw", g=g, neigs=neigs,
                    basisSize=a.compSize, finiteL=True)
            db.insert(datadict=datadict, spec=a.eigenvalues[g]["raw"])

            datadict = dict(k=k, ET=ET, L=L, ren="renloc", g=g, eps=epsraw[g],
                    neigs=neigs, basisSize=a.compSize, finiteL=True)
            db.insert(datadict=datadict, spec=a.eigenvalues[g]["renloc"])

main()
