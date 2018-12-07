from statefuncs3 import *
from sys import argv, exit

m = 1
k = 1

if len(argv) < 3:
    print("{} <Emax> <L>".format(argv[0]))
    exit(1)

Emax = float(argv[1])
L = float(argv[2])

print("Emax={}, L={}".format(Emax,L))

basis = Basis(m, L, k, Emax)
print("Size basis:", len(basis))

for s in basis.stateList:
    print(reprState(s))
