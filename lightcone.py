import statefuncs
import phi4
import renorm
import sys
import scipy
import math
import database
from scipy import array

# saveondb = False
m = 1
# Number of eigenvalues to compute per sector
neigs = 1

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


print("Full basis size: ", a.basis[1].size)


retVEV = []
retmass = []

a.computePotential(k=1)
a.computePotential(k=-1)



# Compute the potential matrices in the low-energy space below ET



for g in glist:

    a.setCouplings(g4=g)
    print("g=", g)

    print("Computing raw eigenvalues")

    k = 1
    a.computeEigval(k, ET, "raw", neigs=neigs)
    v0 = scipy.array(a.eigenvectors["raw"][k][0])
    phi2 = a.V[k][2].M/L
    retVEV.append(phi2.dot(v0).dot(v0))
    print("Norm of v0:", v0.dot(v0))
    print("phi^2 exp value:", retVEV[-1])

    k = -1
    a.computeEigval(k, ET, "raw", neigs=neigs)
    retmass.append(a.spectrum(k=-1,ren="raw")[0])
    print("mass:", retmass[-1])


print(glist)
print(retmass)
print(retVEV)

VEVlist = array(retVEV)
glist = array(glist)

print("Light-cone coupling:")
print([g*(1+12*g*VEV) for g,VEV in zip(glist,VEVlist)])
