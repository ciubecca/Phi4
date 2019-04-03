# Compute eigenvalues for a varying range of ET and g4, for fixed L and g2

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
# Add non-local counterterms
nonlocct = True
# Save lowest eigenvector
eigv = False
# Normalization such that g is divided by 4!
fourfacnorm = False

toy = False

print("ct = {}, eigv={}".format(ct, eigv))

if savedb:
    db = Database()

m = 1
neigs = 4


g4list = array([5*i for i in range(1,21)])/24.
print("g4: ", g4list)

ETmin = 10


# Number of ET's
if toy:
    nET = 1
else:
    nET = 10


print("nET={}".format(nET))

if len(argv) < 4:
    print("{} <L> <ETmax> <g2>".format(argv[0]))
    exit(1)


L = float(argv[1])
Emax = float(argv[2])
g2 = float(argv[3])

print("L={}, ETmax={}, g2={}".format(L, Emax, g2))

ETlist = np.linspace(ETmin, Emax, nET)
print("ETlist: ", ETlist)

t0 = time()
print("Computing basis...")
a = Phi4(m, L, Emax, momcut=False)
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

            a.setmatrix(k, ET)

            print("k={}, ET={}, g2={}".format(k, ET, g2))

            for g4 in g4list:

                if fourfacnorm:
                    a.setg(0, g2, g4/(factorial(4)), ct=ct, cutoff=ET, impr=False)
                else:
                    a.setg(0, g2, g4, ct=ct, cutoff=ET, impr=False)

                # print("Diagonalizing matrix...")
                a.computeEigval(k=k, neigs=neigs, eigv=eigv, nonlocct=nonlocct)
                # print("Spectrum: ", a.eigval)

                print("k={}, {}".format(k,a.eigval[k]))

                if savedb:
                    data = {"neigs":neigs, "logct":ct, "g2":g2, "g4":g4,
                            "spec":a.eigval[k], "L":L, "ET":ET, "Lambda":np.inf, "m":m, "k":k,
                            "momcut":False, "impr":False, "fourfacnorm":fourfacnorm}
                    if eigv:
                         # Save the lowest eigenvector
                        data['eigv'] = a.eigv[k][:, 0]
                    db.insert(data)
    t1 = time()
    print("Elapsed: ",t1-t0)
