import inspect
import os
from phi4 import Phi4
import sys
import math
import scipy
from statefuncs import *


k = +1
m = 1
g = 1


argv = sys.argv

args = "<L> <ET> <row>"
if len(argv) < 4:
    print("python", argv[0], args)
    sys.exit(-1)

L = float(argv[1])
ET = float(argv[2])
row = int(argv[3])

a = Phi4(m,L)
a.buildBasis(Emax=ET)

# vectorlist = a.basis[k].stateList[:10]
vectorlist = [[],
        [(0,1),(1,4),(-4,1)]
        ]
basisl = Basis(k, vectorlist, a.basis[k].helper)
a.basisl[k] = basisl


a.computePotential(k)
a.computeLEVs(k)

V = a.Vll[k][6].M.todense()

print(V)
