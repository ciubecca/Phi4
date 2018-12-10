from phi4 import *
from sys import argv, exit
from math import factorial
from statefuncs import *

m = 1
neigs = 6

if len(argv) < 4:
    print("{} <Emax> <L> <g4> [<g2>]".format(argv[0]))
    exit(1)

Emax = float(argv[1])
L = float(argv[2])
g4 = float(argv[3])
try:
    g2 = float(argv[4])
except IndexError:
    g2 = 0


bases = Basis.fromScratch(m, L, Emax)

for k in (-1,1):
    print("k={}, Emax={}, L={}, size={}".format(k, Emax, L, len(bases[k])))

eigs = {}

for k in (-1,1):
    a = Phi4(bases[k])
    a.computePotential()

    a.setg(0, g2, g4/(factorial(4)))

    a.computeEigval(neigs=neigs)
    eigs[k] = a.eigval

    if k==1:
        # print("state[15]: ", toCanonical(a.basis.stateList[15]))
        # print("state[3]: ", toCanonical(a.basis.stateList[3]))
        # print(a.V[4])
        pass

print("Emax={}, L={}, g2={}, g4={}".format(Emax, L, g2, g4))
print("Vacuum: ", eigs[1][0])
print("Even spectrum: ", eigs[1][1:]-eigs[1][0])
print("Odd  spectrum: ", eigs[-1]-eigs[1][0])
