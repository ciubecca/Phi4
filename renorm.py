import sys
import math
from math import log, pi
from scipy import integrate
import scipy

def ft0(g2, g4, E, m):
    return (g2**2./pi + g4**2.*(-3/(2*pi) + 18/pi**3 * log(E/m)**2))/(E**2.)


def ft2(g2, g4, E, m):
    return (g2*g4*12./pi + g4**2.*72./pi**2. * log(E/m))/(E**2.)


def ft4(g2, g4, E, m):
    return (g4**2. * 36./pi) / (E**2)


def renlocal(g4, EL, m, eps, g2=0):
    deltag = {}

    deltag[2] = {}
    deltag[2][0] = integrate.quad(lambda E: ft0(g2,g4,E,m)/(eps-E),EL,scipy.inf)[0]
    deltag[2][2] = integrate.quad(lambda E: ft2(g2,g4,E,m)/(eps-E),EL,scipy.inf)[0]
    deltag[2][4] = integrate.quad(lambda E: ft4(g2,g4,E,m)/(eps-E),EL,scipy.inf)[0]

    return deltag
