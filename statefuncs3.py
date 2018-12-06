import bisect
from sys import getsizeof as sizeof
import scipy
from scipy import array, pi, sqrt
from math import floor, factorial
from collections import Counter
from itertools import combinations
import itertools
import numpy as np

class Helper():
    """ This is just a "helper" class used to conveniently compute energies of
    oscillators and states and so on"""

    def __init__(self, m, L):
        self.L = L
        self.m = m

    def energy(self, state):
        """ Computes energy of state in Repr1 """
        return sum(Zn*self.omega(n) for n,Zn in state)

    def totwn(self, state):
        return np.sum([Zn*n for n,Zn in state])

    def omega(self, n):
        m = self.m
        return sqrt(m**2 + self.kSq(n))

    def kSq(self, n):
        """ Squared momentum corresponding to wave number n """
        L = self.L
        return (2*pi/L)**2*np.sum(n**2)


class Basis():
    """ Class used to store and compute a basis of states"""

    def __init__(self, m, L, k, Emax, Lambda=np.inf):
        """ Builds the truncated Hilbert space up to cutoff Emax from scratch, in repr1
        m: mass
        L: side of the torus
        Emax: maximal energy of the states
        """

        self.helper = Helper(m, L)
        helper = self.helper
        energy = helper.energy
        self.Emax = Emax
        m = helper.m

        self.NEwnlist = self.genNEwnlist(Emax, Lambda)
        self.stateList = self.buildBasis()
        self.elist = [energy(s) for s in self.stateList]


    def reprState(self, state):
        return [(tuple(n), Zn) for n,Zn in state]

    def __repr__(self):
        return str([self.reprState(s) for s in self.stateList])

    def genNEwnlist(self, Emax, Lambda):
        """ Generate list of North-East moving wave numbers momenta, nx > ny >= 0,
        sorted in energy, below cutoffs Emax and Lambda """

        helper = self.helper
        omega = helper.omega
        kSq = helper.kSq

        ret = []

        for ny in itertools.count():
            for nx in itertools.count(1):

                n = array([nx,ny])

                if omega(n)>Emax/2 or kSq(n)>Lambda**2:
                    if nx==1:
                        ret.sort(key=lambda n: omega(n))
                        return ret
                    else:
                        break

                ret.append(n)

        raise RuntimeError("Shouldn't get here")


    def genNEstatelist(self, NEstate=[], idx=0):
        """ Recursive function generating all North-East moving states in Repr 1 starting from NEstate, by adding
        any number of particles with momentum self.NEwnlist[idx] """

        helper = self.helper
        m = helper.m
        Emax = self.Emax
        omega = helper.omega

        if idx == len(self.NEwnlist):
            return [NEstate]

        # Two-dimensional NE-moving wave number
        n = self.NEwnlist[idx]

        # Keep track of current energy
        E = helper.energy(NEstate)
        # Keep track of current total wave number
        WN = helper.totwn(NEstate)

        ret = []

        for Zn in itertools.count():

            newstate = NEstate[:]

            if Zn>0:
                mode = (n, Zn)

                E += omega(n)
                WN += n

                # We need to add another particle to have 0 total momentum.
                if E + omega(WN) > Emax:
                    break

                newstate.append(mode)

            ret += self.genNEstatelist(newstate, idx+1)

        return ret



    def buildBasis(self):
        """ Generates the basis starting from the list of RM states, in repr1 """

        helper = self.helper
        omega = helper.omega
        m = helper.m
        energy = helper.energy


        NEstatelist = self.genNEstatelist()
        statelist = sorted(NEstatelist, key=energy)

        return statelist
