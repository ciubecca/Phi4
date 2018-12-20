import scipy
from scipy import array, pi, sqrt
from math import floor, factorial
import numpy as np
from symmetry import *

tol = 10**-10

def sortOsc(s):
    """ Sort modes in a state according to momenta """
    return list(sorted(s, key=lambda x: tuple(x[0])))

def occn(s):
    """ Occupation number of state """
    return sum([Zn for n,Zn in s])

# XXX Is this necessary?
def toCanonical(state):
    """ Transorm the state in representation 1 to canonical ordering of the momenta """
    return list(sorted(((tuple(n), Zn) for n,Zn in state), key=lambda x: x[0]))


class Helper():
    """ This is just a "helper" class used to conveniently compute energies of
    oscillators and states and so on"""

    def __init__(self, m, L, Emax, Lambda=np.inf, noscmax=8):
        """ noscmax: max number of oscillators """

        self.L = L
        self.m = m
        self.Emax = Emax
        self.Lambda = Lambda

        # Maximum of sqrt(nx^2 + ny^2)
        self.nmaxFloat = L/(2*pi)*min(Lambda, sqrt((Emax/2)**2-m**2))+tol
        # Maximum integer wave number
        nmax = floor(self.nmaxFloat)
        self.nmax = nmax

        # Do not shift indices, but use negative indices for slight optimization, which avoids some operations in omega()
        self.omegaMat = np.zeros(shape=(2*nmax+1,2*nmax+1))
        for nx in range(-nmax,nmax+1):
            for ny in range(-nmax,nmax+1):
                self.omegaMat[nx][ny] = self._omega(array([nx,ny]))

        # Dictionary of allowed momenta ((kx,ky) -> idx), where idx is used as index of states in representation 2
        self.allowedWn = dict()
        idx = 0
        for nx in range(-nmax, nmax+1):
            for ny in range(-nmax, nmax+1):
                if sqrt(nx**2+ny**2) <= self.nmaxFloat+tol:
                    self.allowedWn[(nx,ny)] = idx
                    idx += 1


        # Set of allowed momenta in first and second quadrants, plus zero momentum
        self.allowedWn12 = set()
        for nx in range(-nmax, nmax+1):
            if nx < 0:
                nymin = 1
            else:
                nymin = 0
            for ny in range(nymin, nmax+1):
                if sqrt(nx**2+ny**2) <= self.nmaxFloat+tol:
                    self.allowedWn12.add((nx,ny))

        occmax = floor(Emax/m+tol)
        self.normFactors = scipy.zeros(shape=(noscmax+1,noscmax+1,occmax+1))
        for c in range(noscmax+1):
            for d in range(noscmax+1):
                for n in range(occmax+1):
                    if d <= n:
                        self.normFactors[c,d,n] = \
                            sqrt(factorial(n)*factorial(n-d+c)/factorial(n-d)**2)
                    else:
                        self.normFactors[c,d,n] = scipy.nan

        self.transfMat = self._genTransfMatrix()

    def torepr2(self, s):
        ret = [0]*len(self.allowedWn)
        for n,Zn in s:
            ret[self.allowedWn[tuple(n)]] = Zn
        return ret

    def oscEnergy(self, wnlist):
        """ Energy of an oscillator (ordered tuple of momenta) """
        return sum(self.omega(n) for n in wnlist)

    # This is slow
    def energy(self, state):
        """ Computes energy of state in Repr1 """
        return sum(Zn*self.omega(n) for n,Zn in state)

    def totwn(self, state):
        if state==[]:
            return array([0,0])
        return sum([Zn*np.array(n) for n,Zn in state])

    def _omega(self, n):
        """ Energy corresponding to wavenumber n"""
        return sqrt(self.m**2 + self.kSq(n))

    def omega(self, n):
        """ Energy corresponding to wavenumber n"""
        return self.omegaMat[n[0]][n[1]]

    # XXX This function can be improved. Need to find the configuration of allowed momenta with minimal energy
    def minEnergy(self, wn, z0min=0):
        """ This function provides a lower bound on the minimal energy to add to a state with total wavenumber "wn",
        in order to create a state with total zero momentum
        z0min: minimum number of particles that must be added """

        m = self.m

        if wn[0]==0 and wn[1]==0:
            return m*z0min
        else:
            # XXX This cannot be used if Lambda < inf
            # return self.omega(wn)
            return self._omega(wn)

    def kSq(self, n):
        """ Squared momentum corresponding to wave number n """
        L = self.L
        return (2*pi/L)**2*(n[0]**2+n[1]**2)

    def maxmom(self, s):
        """ Occupation number of state """
        if s == []:
            return 0.
        return max(sqrt(self.kSq(n)) for n,_ in s)

# TODO Inspect: is it more efficient to transform states to Repr 1 and back, instead?
    def _genTransfMatrix(self):
        """ Generate transformation matrix for all 8 symmetry elements
        for states in representation 2 """

        allowedWn = self.allowedWn
        l = len(allowedWn)
        ret = []

        for op in ((Id, rot, rot2, rot3, refly, reflx, xs, ys)):
            mat = np.zeros(shape=(l,l))

            for wn,i in allowedWn.items():
                # wn = array(wn)
                j = allowedWn[tuple(np.dot(op,array(wn)))]
                mat[i, j] = 1

            ret.append(mat)
        return ret


    def genMomentaPairs(self, totpairsmomenta=None):
        """ Generate sets of all inequivalent pairs of momenta,
        ordered lexicographically, and indexed by total momentum.
        The list of total momenta is
        This is a subroutine used to construct the V22 matrix.
        ksqmax: if we generate oscillators between two bases b1, b2 such that
        E1 < E2, not all pairs of momenta for the annihilation operators are allowed.
            """

        omega = self.omega
        minEnergy = self.minEnergy
        allowedWn = self.allowedWn
        Emax = self.Emax

        # Sort 2d momenta lexicographically
        allowedWnList = list(map(lambda x:np.array(x), sorted(allowedWn)))

        l = len(allowedWnList)
        elist = [omega(wn) for wn in allowedWnList]

        allowedWnPairs = {}

        for i1 in range(l):
            k1 = allowedWnList[i1]
            e1 = elist[i1]

            for i2 in range(i1,l):
                k2 = allowedWnList[i2]
                k12 = k1+k2

                k12 = tuple(k12)

                # ksqmax is the maximum total momentum of pairs of momenta
                # that can be annihilated
                if totpairsmomenta!=None and k12 not in totpairsmomenta:
                    continue

                e12 = e1+elist[i2]

                # XXX CHECK if I can comment this
                # if k12 not in allowedWn:
                    # continue

                if e12+minEnergy(k12) > Emax+tol:
                    continue

                if k12 not in allowedWnPairs:
                    allowedWnPairs[k12] = []

                allowedWnPairs[k12].append((tuple(k1),tuple(k2)))

# Sort lexicographically
        # for k12 in allowedWnPairs.keys():
            # allowedWnPairs[k12] = list(sorted(allowedWnPairs[k12]))

        return allowedWnPairs
