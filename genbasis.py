from phi4 import *
from sys import argv, exit
from math import factorial
from statefuncs import *
import copy
from scipy import sparse
from database import *
from time import time


m = 1

if len(argv) < 5:
    print("{} <L> <Emax> <Lambda> <occmax>".format(argv[0]))
    exit(1)


L = float(argv[1])
Emax = float(argv[2])
Lambda = float(argv[3])
occmax = int(argv[4])

t0 = time()
print("Computing basis...")
bases = Basis.fromScratch(m, L, Emax, Lambda, occmax=occmax)
t1 = time()
print("Elapsed: ",t1-t0)

for k in (-1,1):
    print("k={}, L={}, Emax={}, Lambda={}, occmax={}, size={}".format(k, L, Emax, Lambda, occmax, len(bases[k])))
