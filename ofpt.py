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

print("L={}, ET={}, Lambda={}".format(L,ET,Lambda))

t0 = time()
print("Computing basis...")
bases = Basis.fromScratch(m, L, ET, Lambda)
t1 = time()
print("Elapsed: ",t1-t0)


Vlist = None
V22 = None

subidx = [0]

EL = 2*ET
ELp = EL

basis = bases[1]

print("Generating HE basis")
genHEBasis(basis, subidx, EL, ELp)

for k in (-1,1):
    t0 = time()
    print("Computing k={} matrices...".format(k))
    a = Phi4(bases[k])
    Vlist, V22 = a.computePotential(Vlist, V22)
    t1 = time()
    print("Elapsed: ",t1-t0)
