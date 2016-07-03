from phi1234 import Phi1234
import math
from math import factorial
import scipy


def findEmax(mu, L, maxdim, start, step=1.):
    Emax = start

    a = Phi1234()
    a.buildFullBasis(k=1, L=L, m=mu, Emax=Emax)
    size = a.fullBasis[1].size
    
    while(a.fullBasis[1].size < maxdim):
        size = a.fullBasis[1].size
        Emax += step
        a = Phi1234()
        a.buildFullBasis(k=1, L=L, m=mu, Emax=Emax)

    Emax -= step
    return Emax
