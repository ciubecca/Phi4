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


bases = Basis.fromScratch(m, L, ET, Lambda)
k = 1
maxmom = bases[k].helper.maxmom

# Total overlap with zero modes
zovlp = []
for j,g4 in enumerate(g4list):
    zovlp.append(sum([eigv[k][j][i]**2 for i,s in enumerate(bases[k].stateList) if
        abs(maxmom(s))<tol]))

print(zovlp)

exit(0)

# Mass
plt.figure(1, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
plt.title(title)
plt.xlabel(r"$E_T$")
plt.ylabel(r"$m_{\rm ph}$")
plt.legend(loc=loc)
plt.savefig("plots/massvsET_{}.{}".format(fname,form))
plt.clf()
