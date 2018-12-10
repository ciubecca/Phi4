from phi4 import *
from sys import argv, exit

m = 1

if len(argv) < 4:
    print("{} <Emax> <L> <g2>".format(argv[0]))
    exit(1)

Emax = float(argv[1])
L = float(argv[2])
g2 = float(argv[3])

bases = Basis.fromScratch(m, L, Emax)

for k in (-1,1):
    print("k={}, Emax={}, L={}, size={}".format(k, Emax, L, len(bases[k])))

eigs = {}

for k in (-1,1):
    a = Phi4(bases[k])
    a.computePotential()

    a.setg(0, g2, 0)

    a.computeEigval()
    eigs[k] = a.eigval

# This should be proportional to g2
print(((eigs[-1][0] - eigs[1][0])**2-1)/2)
