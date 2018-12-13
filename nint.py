import math
import scipy
import numpy
import vegas
import time
import os
from numpy import prod
from numpy import arctan
from scipy.special import factorial
from integrator import *

nitn=10
neval=5000

print("Computing O(VV) vacuum diagram...")

lamList = [10, 20, 43, 63]

integ = Phi0_1(nitn, neval)

for lam in lamList:
    print("Integrating for lambda={}".format(lam))
    res = integ.do(lam)

    print("Result: {}".format(res))
