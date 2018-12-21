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
    print("{} <L> <ET> <Lambda>".format(argv[0]))
    exit(1)


L = float(argv[1])
ET = float(argv[2])
Lambda = float(argv[3])


EL = 3*ET
ELp = EL

print("L={}, ET={}, Lambda={}, EL={}".format(L,ET,Lambda, EL))

t0 = time()
print("Computing basis...")
bases = Basis.fromScratch(m, L, ET, Lambda)
bases2 = Basis.fromScratch(m, L, EL, Lambda)


subidx = [0]
Vlist = None
V22 = None

E0 = {1:0, -1:m}

for k in (-1,1):
    basisH = genHEBasis(bases[k], subidx, EL, ELp)

    if k==1:
        idxlist = [i for i,s in enumerate(bases2[k].stateList) if occn(s)==4]
    else:
        helper = bases2[k].helper
        idxlist = [i for i,s in enumerate(bases2[k].stateList) if occn(s)==3 or
                (occn(s)==5 and helper.torepr2(s)[helper.allowedWn[(0,0)]]>0)]

    assert basisH.size == len(idxlist)

    V = genVHl(bases[k], subidx, basisH, L)

    a = Phi4(bases2[k])
    Vlist, V22 = a.computePotential(Vlist, V22)

    V2 = subcolumns(subrows(a.V[4], subidx), idxlist)

    prop = 1/(E0[k]-array(basisH.energyList))

    # deltaE = np.dot(np.dot(V2.transpose(), prop), V2)
    deltaE = np.einsum("ij,j,kj", V2.todense(), prop, V2.todense())[0][0]
    deltaE2 = np.einsum("ij,j,kj", V.todense(), prop, V.todense())[0][0]

    assert abs(deltaE-deltaE2)<tol
