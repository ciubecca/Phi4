import sys
import numpy as np
import scipy
import math
from scipy import pi, log, log10, array, sqrt, stats
import database
from sys import exit, argv
from time import time
from paramplots import *

form  = "pdf"

ct = True

klist = (1,-1)
neigs = 4

db = database.Database()

if len(argv) < 6:
    print("{} <L> <ET> <Lambda> <g4> <g2>".format(argv[0]))
    sys.exit(-1)

L = float(argv[1])
ET = float(argv[2])
Lambda  = float(argv[3])
g4  = float(argv[4])
g2  = float(argv[5])

spectrum = {}
eigv = {}

for k in klist:

    approxQuery = {"g4":g4, "g2":g2, "L":L, "ET":ET, "Lambda":Lambda}
    exactQuery = {"k": k, "neigs":neigs, "logct":ct}
    boundQuery = {}

    try:
        spectrum[k] = array(db.getObjList('spec', exactQuery, approxQuery, orderBy="date")[0])
        eigv[k] =  db.getObjList('eigv', exactQuery, approxQuery, orderBy="date")[0]

    except IndexError:
        print("Not found:", exactQuery, approxQuery)
        exit(-1)

print(eigv[1][:10])

exit(0)

# Mass
plt.figure(1, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
plt.title(title)
plt.xlabel(r"$E_T$")
plt.ylabel(r"$m_{\rm ph}$")
plt.legend(loc=loc)
plt.savefig("plots/massvsET_{}.{}".format(fname,form))
plt.clf()
