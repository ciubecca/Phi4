import gc
import statefuncs
import phi4
import renorm
import sys
import scipy
import math
import database
from profile_support import *

loc3 = True

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

EL = ratioELET*ET
ELp = ratioELpET*ET
ELpp = ratioELppELp*ELp

print("EL, ELp, ELpp:", EL, ELp, ELpp)

def main():

    a = phi4.Phi4(m, L)
    a.buildBasis(Emax=ET)

    for k in klist:

        # Compute the potential matrices in the low-energy space below ET
        a.computePotential(k)

        print("k=", k)
        print("Full basis size: ", a.basis[k].size)

        print("Memory before computeLEVs() ", memory_usage())

        basisl = a.basis[k]
        a.computeLEVs(k=k, basisl=basisl, loc3=loc3)

        print("Memory after computeLEVs() ", memory_usage())

        print("Generating high energy basis...")
        a.genHEBasis(k, EL=EL, ELp=ELp, ELpp=ELpp)
        print("Size of HE basis:", a.basisH[k].size)

        print("Memory: ", memory_usage())

main()
