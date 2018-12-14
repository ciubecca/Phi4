from profile_support import *
from phi4 import *
from sys import argv, exit

m = 1

if len(argv) < 4:
    print("{} <L> <Emax> <Lambda>".format(argv[0]))
    exit(1)

L = float(argv[1])
Emax = float(argv[2])
Lambda = float(argv[3])

bases = Basis.fromScratch(m, L, Emax, Lambda)

for k in (-1,1):
    print("k={}, L={}, Emax={}, Lambda={}, size={}".format(k, L, Emax, Lambda, len(bases[k])))


Vlist = None
V22 = None

for k in (-1,1):
    a = Phi4(bases[k])
    a.computePotential(Vlist, V22)
