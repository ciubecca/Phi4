from phi4 import *
from sys import argv, exit
from math import factorial
from statefuncs import *
import copy
from scipy import sparse
from database import *
from time import time

# Save results on database
savedb = True
# Add counterterms
ct = True
# Save lowest eigenvector
eigv = False

print("ct = {}, eigv={}".format(ct, eigv))

if savedb:
    db = Database()

m = 1
neigs = 4

g4list = np.linspace(2,30,15)
print("g4 ;", g4list)

lammin = 4
ETmin = 10

# Number of Lambda's
nlam = 2
# Number of ET's
nET = 2

print("nET={}, nlam={}".format(nET, nlam))

if len(argv) < 5:
    print("{} <L> <ETmax> <Lambdamax> <g2>".format(argv[0]))
    exit(1)


L = float(argv[1])
Emax = float(argv[2])
Lambda = float(argv[3])
g2 = float(argv[4])

print("L={}, ETmax={}, Lambdamax={}, g2={}".format(L, Emax, Lambda, g2))

lamlist = np.linspace(lammin, Lambda, nlam)
ETlist = np.linspace(ETmin, Emax, nET)

t0 = time()
print("Computing basis...")
a = Phi4(m, L, Emax, Lambda, momcut=True)
t1 = time()
print("Elapsed: ",t1-t0)

for k in (-1,1):
    print("k={}, size={}".format(k, len(a.bases[k])))


print("Computing matrices...".format(k))
t0 = time()
a.computePotential()
t1 = time()
print("Elapsed: ",t1-t0)

eigs = {}

for k in (-1,1):

    for ET in ETlist:
        for lam in lamlist:

            a.setmatrix(k, ET, lam)

            print("k={}, ET={}, lam={}, g2={}".format(k, ET, lam, g2))

            for g4 in g4list:

                # FIXME Error: I should change the renormalization constant according to Lambda!
                a.setg(0, g2, g4/(factorial(4)), ct=ct, cutoff=lam, impr=False)

                # print("Diagonalizing matrix...")
                a.computeEigval(k=k, neigs=neigs, eigv=eigv)
                # print("Spectrum: ", a.eigval)

                if savedb:
                    data = {"neigs":neigs, "logct":ct, "g2":g2, "g4":g4,
                            "spec":a.eigval[k], "L":L, "ET":ET, "Lambda":lam, "m":m, "k":k,
                            "momcut":True, "impr":False}
                    if eigv:
                         # Save the lowest eigenvector
                        data['eigv'] = a.eigv[k][:, 0]
                    db.insert(data)
    t1 = time()
    print("Elapsed: ",t1-t0)
