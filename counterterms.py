# File containing functions to compute counterterms

import scipy
import math
from math import factorial
from scipy import exp, pi, array, e, sqrt, log
import numpy as np
from phi4 import *

m = 1.

def ct0ET(ET, En, m):
    """ Full non-local g^2 correction to the vacuum energy density, from the SW formalism """
    if ET-En < 4*m:
        return 0
    else:
        return -24**2*1/(96*(4*pi)**3)*((ET-En)-8*log((ET-En)/4)-16/(ET-En))

def ct0ETnonloc(ET, En, m):
    """ Partial non-local g^2 correction to the vacuum energy density, subtracting full local g^2 energy density """
    return ct0ET(ET, En, m) - ct0ET(ET, 0, m)


def ct2ET(ET, m):
    """ Full g^2 mass correction, computed via phase space integral """
    return -((24)**2)*1/(12*(4*pi)**2)*\
            (3*log(ET/m)+2*log(ET/m-1)-3*log(ET/m-2)-log(64))*0.5


def ct0ET3(ET, m):
    """ Approximation of the infinite volume g^3 vacuum diagram, in the range of small ET """
    return (24**3)*(-2.81243*10**(-7)+3.47232*10**(-8)*ET)


def ct2Lam(Lambda, m):
    """ Full g^2 mass correction with momentum cutoff in infinite volume, computed via Monte Carlo """
    a = 1/(12*(4*pi)**2)
    b = 3.736124473715983
    c = -a * 1.5848415795962967
    return (24**2)*(-a*log(Lambda/(b*m)) + c*m/Lambda)



class exactct():
    """ Class to compute exact counterterms in finite volume and at finite ET """

    def __init__(self, L, ETmax):
        """
        L: side of the torus
        ETmax: maximum ET that will be ever used.
        """
        self.L = L
        self.ETmax = ETmax

        subidx = {k:[0] for k in (1,)}
        E0 = {1:0}

        a = Phi4(m, L, ETmax)
        basis = a.bases[1]
        self.basisH = genHEBases(a.bases, subidx, ETmax, ETmax, k=1)[1]

        self.VlH = genVHl(basis, subidx[1], self.basisH, L)
        self.prop = 1/(-np.array(self.basisH.energyList))

        self.Vhh = genVhh(self.basisH, L)


    def ct2(self, ET, En):
        """ Return the exact non-local counterterm in finite volume
        ET: Energy cutoff
        En: energy of the external state
        """

        if ET-En<4:
            return 0.

        idxlist = self.basisH.subidxlist(ET-En, Emin=2)
        Vsub = subcolumns(self.VlH, idxlist)
        propsub = self.prop[idxlist]
        proj = scipy.sparse.spdiags(propsub, 0, len(propsub), len(propsub)).tocsc()

        deltaE = (Vsub*proj*Vsub.transpose()).todense()[0,0]

        return deltaE/self.L**2

    def ct3(self, ET):
        """ Return exact g^3 local vacuum counterterm in finite volume """

        idxlist = self.basisH.subidxlist(ET, Emin=2)
        VlHsub = subcolumns(self.VlH, idxlist)
        Vhhsub = submatrix(self.Vhh, idxlist)

        propsub = self.prop[idxlist]
        proj = scipy.sparse.spdiags(propsub, 0, len(propsub), len(propsub)).tocsc()

        deltaE = (VlHsub*proj*Vhhsub*proj*VlHsub.transpose()).todense()[0,0]
        return deltaE/self.L**2
