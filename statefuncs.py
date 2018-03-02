import bisect
from sys import getsizeof as sizeof
import scipy
from scipy import array, pi, sqrt
from math import floor, factorial
from collections import Counter
from itertools import combinations
import itertools
import numpy as np

tol = 10**-10

class Helper():
    """ This is just a "helper" class used to conveniently compute energies of
    oscillators and states and so on"""
    def __init__(self,m,L,Emax,noscmax=8):
        self.L = L
        self.m = m
        self.nmax = self.Emaxtonmax(Emax)
        self.wnList = array(range(-self.nmax,self.nmax+1))
        occmax = int(floor(Emax/m))

        self.normFactors = scipy.zeros(shape=(noscmax+1,noscmax+1,occmax+1))
        for c in range(noscmax+1):
            for d in range(noscmax+1):
                for n in range(occmax+1):
                    if d <= n:
                        self.normFactors[c,d,n] = \
                            sqrt(factorial(n)*factorial(n-d+c)/factorial(n-d)**2)
                    else:
                        self.normFactors[c,d,n] = scipy.nan

        self.omegaList = self._omega(self.wnList)


    def energy(self, state):
        """ Computes energy of state in Repr1 """
        return sum(Zn*self.omega(n) for n,Zn in state)

    def energy2(self, state):
        """ Computes energy of state in Repr2 """
        return sum(Zn*self.omega(n-self.nmax) for n,Zn in enumerate(state))

    def totalWN(self, state):
        return sum(Zn*n for n,Zn in state)

    def oscEnergy(self, wnlist):
        """ Energy of an oscillator (ordered tuple of momenta) """
        return sum(self.omega(n) for n in wnlist)



    def calcOscEnergyDict(self):
        clists = []
        nmax = self.nmax

        clists.append(())

        for k1 in range(-nmax,nmax+1):
            clists.append((k1,))

            for k2 in range(k1,nmax+1):
                clists.append((k1,k2))

                for k3 in range(k2,nmax+1):
                    clists.append((k1,k2,k3))

                    for k4 in range(k3,nmax+1):
                        clists.append((k1,k2,k3,k4))

        self.oscEnergyDict = {clist: self.oscEnergy(clist) for clist in clists}

    def _omega(self, n):
        return sqrt(self.m**2+((2*pi/self.L)*n)**2)

    def omega(self, n):
        """ Energy corresponding to wavenumber n"""
        return self.omegaList[n+self.nmax]

    def torepr2(self, state):
        """ Transform state from repr1 to repr2 """
        ret = [0]*(2*self.nmax+1)
        for n,Zn in state:
            ret[n+self.nmax] = Zn
        return ret

    def torepr1(self, state):
        """ Transform state from repr2 to repr1 """
        return [(self.wnList[i],state[i]) for i in range(2*self.nmax+1)
            if state[i]!= 0]

    def Emaxtonmax(self, Emax):
        """ return nmax corresponding to given Emax """
        return int(floor(sqrt((Emax/2.)**2.-self.m**2.)*self.L/(2*pi)))


# Mantains the ordering of wavenumbers
def reverse(state):
    """ Apply the spatial parity transformation to a state in representation 1"""
    return [(-n,Zn) for n,Zn in state[::-1]]

def occn(state):
    """ Computes the occupation number of a state"""
    return sum(Zn for n,Zn in state)


def isSorted(x, key):
    return all(key(x[i]) <= key(x[i+1]) for i in range(len(x)-1))

class Basis():
    """ Class used to store and compute a basis of states"""
    def __init__(self, k, stateset, helper, repr1=True, repr1Emax=None):
        """ Standard constructor
        k: parity quantum number
        stateset: set or list of states in representation 1
        helper: Helper object
        orderEnergy: if True order the list of vectors in energy
        calcPos: if True construct the dictionary of the positions of the vectors
        """
        self.k = k
        self.helper = helper
        self.nmax = self.helper.nmax

        if repr1:
            energy = helper.energy
            self.stateList = sorted(stateset, key=energy)
            self.energyList = [energy(state) for state in self.stateList]
            self.occnList = [occn(state) for state in self.stateList]
            self.parityList = [int(state==reverse(state)) for state in self.stateList]
            self.repr2List = [bytes(helper.torepr2(state)) for state in self.stateList]
            self.repr1 = True
        # We don't transform to repr1 all the states, but only those
        # we'll cycle over
        else:
            energy = helper.energy2
            self.repr2List = sorted(stateset, key=energy)
            self.stateList = [helper.torepr1(state) for state in self.repr2List if
                    energy(state)<=repr1Emax+tol]
            self.occnList = [sum(state) for state in self.repr2List]
            self.parityList = [int(state==state[::-1]) for state in self.repr2List]
            self.energyList = [energy(state) for state in self.repr2List]
            self.repr1 = False


        self.size = len(self.energyList)
        self.Emax = max(self.energyList)
        self.Emin = min(self.energyList)
        self.occmin = min(self.occnList)

    def irange(self, Erange):
        """ Return the min and max indices for states with energy between
        Emin and Emax """
        Emin = Erange[0]
        Emax = Erange[1]+tol
        imin = bisect.bisect_left(self.energyList, Emin)
        imax = bisect.bisect_left(self.energyList, Emax)
        return range(imin, imax)

    def propagator(self, eps, Emin, Emax):
        """ Return the propagator for states between Emin and Emax """
        v = [1/(eps-e) if Emin<e<=Emax else 0 for e in self.energyList]
        # print([e for e in self.energyList  if Emin<e<=Emax ][0])
        # print(eps)
        # print([1/(eps-e) for e in self.energyList  if Emin<e<=Emax ][0])
        return scipy.sparse.spdiags(v, 0, self.size, self.size).tocsc()

    @classmethod
    def fromScratch(self, m, L, k, Emax, occmax=None, occmin=0):
        """ Builds the truncated Hilbert space up to cutoff Emax from scratch
        m: mass
        L: size of the cylinder
        Emax: maximal energy of the states
        occmax: cutoff in occupation number (optional)"""

        self.helper = Helper(m, L, Emax)
        helper = self.helper
        self.Emax = Emax
        m = helper.m

        self.occmax = occmax
        if occmax ==None:
            self._occmax = int(floor(Emax/m))
        else:
            self._occmax = occmax

        self.occmin = occmin

        # self.nmax is the actual maximum occupied wavenumber of the states
        self.nmax = helper.nmax

        bases = self.buildBasis(self)
        return self(k,bases[k],helper)

    def MBsize(self):
        ret = 0
        """ 64 is the size of a tuple of two ints """
        ret += sum(64*len(v) for v in self.stateList)/10**6
        # Add the memory of the auxiliary vectors
        ret += self.size*32/10**6
        return ret

    def __repr__(self):
        return str(self.stateList)

    def genRMlist(self, RMstate=[], n=1):
        """ Recursive function generating all the states starting from RMstate, by adding
        any number of particles with wavenumber n.
        It starts from the seed state with 0 particles and wavenumber 1 """

        if n > self.nmax:
            return [RMstate]

        maxN = int(floor(min(
                self._occmax-occn(RMstate),
                (self.nmax-self.helper.totalWN(RMstate))/n,
                (self.Emax-self.helper.energy(RMstate))/self.helper.omega(n))))
        ret = []
        for Zn in range(maxN+1):
            newstate = RMstate[:]
            if Zn>0:
                newstate.append((n,Zn))
            ret += self.genRMlist(self,newstate,n+1)
        return ret


    def buildBasis(self):
        """ Generates the basis starting from the list of RM states """

        omega = self.helper.omega
        m = self.helper.m

        RMlist = self.genRMlist(self)

        # divides the list of RMstates into a list of lists,
        # so that two states in each list have a fixed total RM wavenumber,
        # and sort each sublist in energy
        sortedRMlist = sorted(RMlist, key=self.helper.totalWN)
        dividedRMlist = [sorted(l, key=self.helper.energy) for wn,l in
            itertools.groupby(sortedRMlist,key=self.helper.totalWN)]
        statelist = {1:[], -1:[]}

        for RMwn, RMsublist in enumerate(dividedRMlist):
            for i, RMstate in enumerate(RMsublist):
                ERM = self.helper.energy(RMstate)
                ORM = occn(RMstate)

                # LM part of the state will come from the same sublist.
                # We take the position of LMState to be greater or equal
                # to the position of RMstate
                for LMstate in RMsublist[i:]:
                    # we will just have to reverse it
                    ELM = self.helper.energy(LMstate)
                    OLM = occn(LMstate)
                    deltaE = self.Emax - ERM - ELM
                    deltaOcc = self._occmax - ORM - OLM

                    # if this happens, we can break since subsequent LMstates
                    # have even higher energy (RMsublist is ordered in energy)
                    if deltaE < 0:
                        break
                    if deltaOcc < 0:
                        continue

                    maxN0 = int(floor(min(deltaOcc, deltaE/m)))

                    # possible values for the occupation value at rest
                    for N0 in range(maxN0+1):

                        # Condition on minimum occupation number
                        if ORM+OLM+N0 < self.occmin:
                            continue

                        # Only states with correct parity
                        if N0==0:
                            state = reverse(LMstate)+RMstate
                        else:
                            state = reverse(LMstate)+[(0,N0)]+RMstate

                        statelist[(-1)**(N0+OLM+ORM)].append(state)

        return statelist
