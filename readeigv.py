# Print eigenvectors from database, and make a plot of the zero mode contributions

import sys
import numpy as np
import scipy
import math
from scipy import pi, log, log10, array, sqrt, stats
import database
from sys import exit, argv
from time import time
from paramplots import *
from statefuncs import *

form  = "pdf"

ct = True

klist = (1,-1)
neigs = 4
m = 1

db = database.Database()

g4list = np.linspace(2,30,15)

if len(argv) < 5:
    print("{} <L> <ET> <Lambda> <g2>".format(argv[0]))
    sys.exit(-1)

L = float(argv[1])
ET = float(argv[2])
Lambda  = float(argv[3])
g2  = float(argv[4])

spectrum = {k:[] for k in klist}
eigv = {k:[] for k in klist}

for k in klist:
    for g4 in g4list:

        approxQuery = {"g4":g4, "g2":g2, "L":L, "ET":ET, "Lambda":Lambda}
        exactQuery = {"k": k, "neigs":neigs, "logct":ct}
        boundQuery = {}

        try:
            spectrum[k].append(array(db.getObjList('spec', exactQuery, approxQuery, orderBy="date")[0]))
            eigv[k].append(db.getObjList('eigv', exactQuery, approxQuery, orderBy="date")[0])

        except IndexError:
            print("Not found:", exactQuery, approxQuery)
            exit(-1)


print("Computing basis...")
bases = Basis.fromScratch(m, L, ET, Lambda)
k = 1
maxmom = bases[k].helper.maxmom

# Total overlap with zero modes
zovlp = []
for j,g4 in enumerate(g4list):
    zovlp.append(sum([eigv[k][j][i]**2 for i,s in enumerate(bases[k].stateList) if
        abs(maxmom(s))<tol]))

print(zovlp)

plt.figure(1)
plt.plot(g4list, zovlp)

title = r"$g_2$={}, $L$={}, $E_T$={}, $\Lambda$={}".format(g2, L, ET, Lambda)
fname = "g2={}_L={}_ET={}_Lambda={}".format(g2, L, ET, Lambda)


plt.figure(1, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
plt.title(title)
plt.xlabel(r"$g$")
plt.ylabel(r"$\sum_{n_0} \mid \langle \psi \mid n_0 \rangle \mid^2$")
plt.savefig("zmodeovlp_{}.{}".format(fname,form))
plt.clf()
