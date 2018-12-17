from phi4 import *
from sys import argv, exit
from math import factorial
from statefuncs import *
from scipy import sparse
from database import *
from time import time


L = 10
Lambda = np.inf
occmax = 4

ETmin= 10
ETmax= 20
ETlist = np.linspace(ETmin, Emax, 10)

bases = Basis.fromScratch(m, L, ETmax, Lambda, occmax=occmax)

print("L={}, Emax={}, Lambda={}, occmax={}".format(L, ETmax, Lambda, occmax))

for k in (-1,1):
    print("k={}, size={}".format(k, len(bases[k])))

eigs = {}

Vlist = None
V22 = None

k = 1
a = Phi4(bases[k])
a.computePotential()

ret = []

for ET in ETlist:
    a.setmatrix(Emax=ET)
    el = array([e for e in bases[k].energyList if e<=ET+tol])
    prop = np.eye(1/(-el))
    V = a.Vcomp[4]
    ret.append(np.dot(V,prop).dot(V))
