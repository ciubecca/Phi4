import statefuncs
import phi4
import renorm
import sys
import scipy
import math


m = 1
# List of parity quantum numbers
klist = (1,)

# Ratio between EL and ET
ratioELET = 3
# Ratio between ELp and ET
ratioELpET = 1.5
# Ratio between ELpp and ELp
ratioELppELp = 1.5

neigs = 10

argv = sys.argv
if len(argv) < 5:
    print(argv[0], "<L> <g> <ET> <ntails>")
    sys.exit(-1)

L = float(argv[1])
g = float(argv[2])
ET = float(argv[3])
ntails = int(argv[4])

EL = ratioELET*ET
ELp = ratioELpET*ET
ELpp = ratioELppELp*ELp

print("EL/ET:", ratioELET)
print("ELp/ET:", ratioELpET)
print("ELpp/ELp", ratioELppELp)


a = phi4.Phi4(m, L)
a.buildBasis(Emax=ET)

for k in klist:
    print("k=", k)

    print("Full basis size: ", a.basis[k].size)
    a.computePotential(k)

# Computing the high energy bases dimension for highest Emax
    a.setCouplings(g4=g)
    a.computeEigval(k, ET, "raw")


    # Select a set of tails and construct a Basis object, ordered in overlap with
    # the vacuum
    vectorlist = [state for i,state in sorted(enumerate(a.basis[k]), key=lambda x:
            -abs(a.eigenvectors["raw"][k][0][x[0]]))][:ntails]
    basisl = statefuncs.Basis(k, vectorlist, a.basis[k].helper, orderEnergy=False)
    print("ntails:", basisl.size)


    print("Generating high energy basis for highest Emax...")
    # Generate the high-energy "selected" basis by passing a set of tails
    # and a maximum cutoff EL
    a.genHEBases(k, basisl, EL=EL, ELpp=ELpp)
    print("Size of HE basis for DH2:", a.basisH[k].size)
    print("Size of HE basis for DH3:", a.basish[k].size)
