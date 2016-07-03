import sys
import math
from math import log, pi
from scipy import integrate
import scipy

def ft0(g2, g4, E, m, cutoff=0.):
    if E<cutoff:
        return 0.
    else:
        return (g2**2./pi + g4**2.*(-3/(2*pi) + 18/pi**3 * log(E/m)**2))/(E**2.)


def ft2(g2, g4, E, m, cutoff=0.):
    if E<cutoff:
        return 0.
    else:
        return (g2*g4*12./pi + g4**2.*72./pi**2. * log(E/m))/(E**2.)


def ft4(g2, g4, E, m, cutoff=0.):
    if E<cutoff:
        return 0.
    else:
        return (g4**2. * 36./pi) / (E**2)


def renlocal(g0, g2, g4, Emax, Er, m):
    # m is for finite L
    g0r = g0 - integrate.quad(lambda E: ft0(g2,g4,E,m)/(E-Er),Emax,scipy.inf)[0]
    g2r = g2 - integrate.quad(lambda E: ft2(g2,g4,E,m)/(E-Er),Emax,scipy.inf)[0]
    g4r = g4 - integrate.quad(lambda E: ft4(g2,g4,E,m)/(E-Er),Emax,scipy.inf)[0]

    return (g0r,g2r,g4r)

def rensubl(g2, g4, Ebar, Emax, Er, m, cutoff):
    ktab = scipy.linspace(0.,Emax,30,endpoint=True)
    rentab = [[0.,0.,0.] for k in ktab]

    for i,k in enumerate(ktab):
        g0subl = - integrate.quad(lambda E: ft0(g2,g4,E-k,m,cutoff)/(E-Ebar) - ft0(g2,g4,E,m,cutoff)/(E-Er),Emax,scipy.inf)[0]
        g2subl = - integrate.quad(lambda E: ft2(g2,g4,E-k,m,cutoff)/(E-Ebar) - ft2(g2,g4,E,m,cutoff)/(E-Er),Emax,scipy.inf)[0]
        g4subl = - integrate.quad(lambda E: ft4(g2,g4,E-k,m,cutoff)/(E-Ebar) - ft4(g2,g4,E,m,cutoff)/(E-Er),Emax,scipy.inf)[0]
        
        rentab[i] = [g0subl, g2subl, g4subl]

    return ktab,scipy.array(rentab)

