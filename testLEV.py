import statefuncs
from  profile_support import *
import phi4
import renorm
import sys
import scipy
import math
import database


memdbg = False
if memdbg:
    warnings.warn("Running with memory debugging")

m = 1
neigs = 6
glist = [1]


argv = sys.argv
if len(argv) < 4:
    print(argv[0], "<k> <L> <ET>")
    sys.exit(-1)

k = int(argv[1])
L = float(argv[2])
ET = float(argv[3])


print("k, L, ET", k, L, ET)


@profile
def main():


    a = phi4.Phi4(m, L)
    a.buildBasis(Emax=ET)

    # Compute the potential matrices in the low-energy space below ET
    a.computePotential(k)

    print("k=", k)
    print("Full basis size: ", a.basis[k].size)

    a.setglist(glist=glist)


# Compute the raw eigenvalues for cutoff ET
    a.computeEigval(k, ET, "raw", neigs=neigs)
    epsraw = {g: a.eigenvalues[g]["raw"][k][0] for g in glist}
    print("Raw vacuum:", epsraw)

    # Always consider maximal set of tails
    basisl = a.basis[k]
    print("Total number of tails:", basisl.size)

    if memdbg:
        print("memory taken before computeLEVs", memory_usage())

    a.computeLEVs(k, basisl)

    if memdbg:
        print("memory taken before genHEBasis", memory_usage())


main()
