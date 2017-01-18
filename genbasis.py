import statefuncs
import phi4
import renorm
import sys
import scipy
import math
import database

# saveondb = False
m = 1
# Number of eigenvalues to compute per sector
neigs = 1
# List of parity quantum numbers
klist = (1,)
maxntails = 300

# Ratio between EL and ET
ratioELET = 3
# Ratio between ELp and ET
ratioELpET = 2
# Ratio between ELpp and ELp
ratioELppELp = 2


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


@profile
def main():

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
        vectorlist = [state for i,state in sorted(enumerate(a.basis[k]), key=lambda x:
                -abs(a.eigenvectors["raw"][k][0][x[0]]))][:maxntails]
        basisl = statefuncs.Basis(k, vectorlist, a.basis[k].helper)
        print("Total number of tails:", basisl.size)


        print("Generating high energy basis...")
        # Generate the high-energy "selected" basis by passing a set of tails
        # and a maximum cutoff EL
        a.genHEBasis(k, basisl, EL=EL, ELp=ELp, ELpp=ELpp)
        print("Size of HE basis:", a.basisH[k].size)


main()
