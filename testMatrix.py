import inspect
import os
from phi4 import Phi4
from newphi4Lorenzo import newPhi4
import sys
import math
import scipy
from statefuncs import *


k = +1
m = 1
g = 1


argv = sys.argv

args = "<L> <ET>"
if len(argv) < 3:
    print("python", argv[0], args)
    sys.exit(-1)

L = float(argv[1])
ET = float(argv[2])

a = newPhi4()
a.buildBasisAndMatrix(L,ET,m)
V1 = a.potential[k]

b = Phi4(m,L)
b.buildBasis(Emax=ET)
b.computePotential(k)
V2 = b.V[k][4].M

print(type((V1-V2).tocoo()))
print(max(abs((V1-V2).tocoo().data)))
