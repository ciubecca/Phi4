from statefuncs3 import *
from sys import argv, exit

m = 1
k = 1

if len(argv) < 4:
    print("{} <k> <Emax> <L>".format(argv[0]))
    exit(1)

k = int(argv[1])
Emax = float(argv[2])
L = float(argv[3])

basis = Basis(m, L, k, Emax)

# for s in basis.stateList:
    # print(reprState(s))

# for i,e in enumerate(basis.elist):
    # print("{}: {}".format(i,e))


print("k={}, Emax={}, L={}, size={}".format(k, Emax,L,len(basis)))
