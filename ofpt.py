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
t1 = time()
print("Elapsed: ",t1-t0)


Vlist = None
V22 = None

for k in (-1,1):
    subidx = [0]
    basis1 = genHEBasis(bases[k], subidx, EL, ELp)

    if k==1:
        sl2 = [s for s in bases2[k].stateList if occn(s)==4]
    else:
        helper = bases2[k].helper
        sl2 = [s for s in bases2[k].stateList if occn(s)==3 or
                (occn(s)==5 and helper.torepr2(s)[helper.allowedWn[(0,0)]]>0)]



    assert basis1.size == len(sl2)

exit(0)

for k in (-1,1):
    t0 = time()
    print("Computing k={} matrices...".format(k))
    a = Phi4(bases2[k])
    Vlist, V22 = a.computePotential(Vlist, V22)
    t1 = time()
    print("Elapsed: ",t1-t0)
