import scipy
import math
from math import factorial
from scipy import exp, pi, array, e, sqrt, log
import numpy as np

def ct0ET(ET, En, m):
    """ Full non-local g^2 correction to the vacuum energy density, from the SW formalism """
    if ET-En < 4*m:
        return 0
    else:
# FIXME Should there be a factor 1/m here?
        return -24**2*1/(96*(4*pi)**3)*((ET-En)-8*log((ET-En)/4)-16/(ET-En))

def ct0ETnonloc(ET, En, m):
    """ Partial non-local g^2 correction to the vacuum energy density, subtracting full local g^2 energy density """
# NOTE This has dimension 3, like the vacuum energy density should
    return ct0ET(ET, En, m) - ct0ET(ET, 0, m)

# def ct2(ET):
    # return -((24)**2)*1/(12*(4*pi)**2)*(log(ET/4)-3/4 +3/ET)

def ct2ET(ET, m):
    """ Full g^2 mass correction, computed via phase space integral """
    return -((24)**2)*1/(12*(4*pi)**2)*\
            (3*log(ET/m)+2*log(ET/m-1)-3*log(ET/m-2)-log(64))*0.5


def ct2Lam(Lambda, m):
    " Full g^2 mass correction, computed via Monte Carlo """
    a = 1/(12*(4*pi)**2)
    b = 3.736124473715983
    c = -a * 1.5848415795962967
    return (24**2)*(-a*log(Lambda/(b*m)) + c*m/Lambda)

