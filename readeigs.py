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

if len(argv) < 6:
    print("{} <L> <ET> <Lambda> <g4> <g2>".format(argv[0]))
    sys.exit(-1)

L = float(argv[1])
ET = float(argv[2])
Lambda  = float(argv[3])
g4  = float(argv[4])
g2  = float(argv[5])

spectrum = {}
masses = {}

for k in klist:

        approxQuery = {"g4":g4, "g2":g2, "L":L, "ET":ET, "Lambda":Lambda}
        exactQuery = {"k": k, "neigs":neigs, "logct":ct}
        boundQuery = {}

        try:
            spectrum[k] = array(db.getObjList('spec', exactQuery, approxQuery, orderBy="date")[0])

        except IndexError:
            print("Not found:", exactQuery, approxQuery)
            exit(-1)

masses[1] = spectrum[1][1:]-spectrum[1][0]
masses[-1] = spectrum[-1]-spectrum[1][0]

print(masses[1])
print(masses[-1])
