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

glist = scipy.linspace(0,3,16)

print("glist", glist)

argv = sys.argv
if len(argv) < 3:
    print(argv[0], "<L> <ET>")
    sys.exit(-1)

L = float(argv[1])
ET = float(argv[2])


a = phi4.Phi4(m, L)
a.buildBasis(Emax=ET)


ret = []

for k in klist:

    # Compute the potential matrices in the low-energy space below ET
    a.computePotential(k)

    print("k=", k)

    print("Full basis size: ", a.basis[k].size)

    for g in glist:

        a.setCouplings(g4=g)
        print("g=", g)


        print("Computing raw eigenvalues")
        a.computeEigval(k, ET, "raw", neigs=neigs)


        print("Raw vacuum:", a.eigenvalues["raw"][k][0])

        v0 = scipy.array(a.eigenvectors["raw"][k][0])
        V2 = a.V[k][2].M
        print(V2.shape, v0.shape)
        ret.append(V2.dot(v0).dot(v0))
        print("phi^2 exp value:", ret[-1])


print(glist)
print(ret)
