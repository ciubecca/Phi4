from statefuncs3 import *
from sys import argv, exit

m = 1

if len(argv) < 3:
    print("{} <Emax> <L>".format(argv[0]))
    exit(1)

Emax = float(argv[1])
L = float(argv[2])

bases = Basis.fromScratch(m, L, Emax)

# for s in basis.stateList:
    # print(reprState(s))

# for i,e in enumerate(basis.elist):
    # print("{}: {}".format(i,e))

for k in (-1,1):
    print("k={}, Emax={}, L={}, size={}".format(k, Emax, L, len(bases[k])))
