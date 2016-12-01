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

    def totalWN(self, state):
        return sum(Zn*n for n,Zn in state)

    def oscEnergy(self, wnlist):
        """ Energy of an oscillator (ordered tuple of momenta) """
        return sum(self.omega(n) for n in wnlist)

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


class Basis():
    """ Class used to store and compute a basis of states"""
    def __init__(self, k, stateset, helper, orderEnergy=True):
        """ Standard constructor
        k: parity quantum number
        stateset: set or list of states in representation 1
        helper: Helper object
        """
        self.k = k
        self.helper = helper
        energy = helper.energy

        self.orderEnergy = orderEnergy
        # Order the states in energy
        if orderEnergy:
            self.stateList = list(sorted(stateset, key=energy))
        else:
            self.stateList = list(stateset)
        self.size = len(self.stateList)

        self.energyList = [energy(state) for state in self.stateList]
        self.occnList = [occn(state) for state in self.stateList]
        self.parityList = [int(state==reverse(state)) for state in self.stateList]

        try:
            self.Emax = max(self.energyList)
            self.Emin = min(self.energyList)
        except IndexError:
            self.Emax = None
            self.Emin = None

        self.statePos = {}
        for i,state in enumerate(self.stateList):
            self.statePos[tuple(helper.torepr2(state))] = i
            self.statePos[tuple(helper.torepr2(state)[::-1])] = i


    @classmethod
    def fromScratch(self, m, L, Emax, occmax=None):
        """ Builds the truncated Hilbert space up to cutoff Emax from scratch
        m: mass
        L: size of the cylinder
        Emax: maximal energy of the states
        occmax: cutoff in occupation number (optional)"""

        self.helper = Helper(m, L, Emax)
        helper = self.helper
        self.Emax = Emax
        m = helper.m

        # This can be "None"
        self.occmax = occmax
        if occmax==None:
            # This is used in the basis construction
            self._occmax = int(floor(Emax/m))
        else:
            self._occmax = occmax

        # self.nmax is the actual maximum occupied wavenumber of the states
        self.nmax = helper.nmax

        bases = self.buildBasis(self)
        return {k:self(k,bases[k],helper) for k in (-1,1)}


    def sub(self, filterFun):
        """ Extracts a sub-basis with vectors v such that filterFun(v)=True """
        return Basis(self.k, filter(filterFun, self.stateList), self.helper,
                self.orderEnergy)


    def lookup(self, state):
        """ Lookup function to look up the position of a state in the basis,
        in representation 1"""
        return self.statePos[tuple(self.helper.torepr2(state))]


    def __repr__(self):
        return str(self.stateList)

    def __getitem__(self,index):
        return self.stateList[index]


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
                        # Only states with correct parity
                        if N0==0:
                            state = reverse(LMstate)+RMstate
                        else:
                            state = reverse(LMstate)+[(0,N0)]+RMstate
                        statelist[(-1)**(N0+OLM+ORM)].append(state)

        return statelist
