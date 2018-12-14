from phi4 import *
from sys import argv, exit
from math import factorial
from statefuncs import *
import copy
from scipy import sparse
from database import *

g2 = -0.75
print("g2 = {}".format(g2))


savedb = True

if savedb:
    db = Database()

m = 1
neigs = 4


if len(argv) < 4:
    print("{} <L> <Emax> <Lambda>".format(argv[0]))
    exit(1)

L = float(argv[1])
Emax = float(argv[2])
Lambda = float(argv[3])

lammin = 4
ETmin = 10
nlam = 2
nET = 10

g4list = np.linspace(0,30,15)
lamlist = np.linspace(lammin, Lambda, nlam)
ETlist = np.linspace(ETmin, Emax, nET)

print("Computing basis...")
bases = Basis.fromScratch(m, L, Emax, Lambda)

for k in (-1,1):
    print("k={}, L={}, Emax={}, Lambda={}, size={}".format(k, L, Emax, Lambda, len(bases[k])))

eigs = {}

Vlist = None
V22 = None

for k in (-1,1):

    print("Computing k={} matrices...".format(k))
    a = Phi4(bases[k])
    Vlist, V22 = a.computePotential(Vlist, V22)

    for ET in ETlist:
        for lam in lamlist:

            a.setmatrix(ET, lam)

            for g4 in g4list:
                a.setg(0, g2, g4/(factorial(4)))
                print("k={}, ET={}, lam={}, g4={}, g2={}".format(k, ET, lam, g4, g2))
                print("Diagonalizing matrix...")
                a.computeEigval(neigs=neigs)
                print("Spectrum: ", a.eigval)

                if savedb:
                    data = {"neigs":neigs, "g2":g2, "g4":g4, "spec":a.eigval, "L":L, "ET":ET, "Lambda":lam, "m":m, "k":k}
                    db.insert(data)
