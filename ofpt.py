from nlo import *
from sys import argv, exit
from math import factorial
from statefuncs import *
import copy
from scipy import sparse
from database import *
from time import time
from phi4 import *
m = 1

if len(argv) < 4:
    print("{} <L> <EL> <Lambda>".format(argv[0]))
    exit(1)


L = float(argv[1])
EL = float(argv[2])
Lambda = float(argv[3])

ET = 2
ELp = EL

print("L={}, Lambda={}, EL={}".format(L, Lambda, EL))

print("Computing basis...")
bases = Basis.fromScratch(m, L, ET, Lambda)


subidx = [0]

E0 = {1:0, -1:m}

for k in (-1,1):
    basisH = genHEBasis(bases[k], subidx, EL, ELp)
    print("k={} basis size={}".format(k, basisH.size))

    V = genVHl(bases[k], subidx, basisH, L)

    prop = 1/(E0[k]-array(basisH.energyList))

    deltaE = np.einsum("ij,j,kj", V.todense(), prop, V.todense())[0][0]
    print(deltaE)
