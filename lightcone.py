import statefuncs
from  profile_support import *
import phi4
import renorm
import sys
import scipy
import math
import database
import numpy as np


m = 1
# Number of eigenvalues to compute per sector
neigs = 6
k =1


argv = sys.argv
if len(argv) < 3:
    print(argv[0], " <L> <ET>")
    sys.exit(-1)

L = float(argv[1])
ET = float(argv[2])


print("L, ET", k, L, ET)

glist = np.linspace(0.1, 2, 11)
print("glist", glist)


def main():

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

main()
