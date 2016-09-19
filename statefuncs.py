import scipy
from scipy import array, pi, sqrt
from math import floor, factorial
from collections import Counter
from itertools import combinations
import itertools
import numpy as np


class phi4Info():
    def __init__(self,m,L,Emax,noscmax=4):
        self.L = L
        self.m = m
        self.Emax = Emax
        self.nmax = int(floor(sqrt((Emax/2.)**2.-m**2.)*L/(2*pi)))
        self.occmax = int(floor(Emax/m))
        self.wnList = array(range(-self.nmax,self.nmax+1))
        self.RMwnList = array(range(1,self.nmax+1))

        self.normFactors = scipy.zeros(shape=(noscmax+1,noscmax+1,self.occmax+1))
        for c in range(noscmax+1):
            for d in range(noscmax+1):
                for n in range(self.occmax+1):
                    if d <= n:
                        self.normFactors[c,d,n] = \
                            sqrt(factorial(n)*factorial(n-d+c)/factorial(n-d)**2)
                    else:
                        self.normFactors[c,d,n] = scipy.nan

        # XXX check
        self.parityFactors = scipy.array([[1,sqrt(2)],[1/sqrt(2),1]])

        self.omegaList = self._omega(self.wnList)
        self.omegaRMList = self._omega(self.RMwnList)

    def computeNorm(self, c, d, n):
        return self.normFactors[c,d,n]

    def energy(self, state):
        return scipy.sum(self.omegaList*state)

    def RMenergy(self, state):
        return scipy.sum(self.omegaRMList*state)

    def RMtotalWN(self, state):
        return sum(self.RMwnList*state)

    def oscEnergy(self, wnlist):
        return sum([self.omega(n) for n in wnlist])

    def _omega(self, n):
        return sqrt(self.m**2+((2*pi/self.L)*n)**2)

    def omega(self, n):
        return self.omegaList[n+self.nmax]

    def repr1torepr2(self, state):
        """ Transform state from repr1 to repr2 """
        return array([state.get(n,0) for n in self.wnList])

    def repr2torepr1(self, state):
        return Counter({self.wnList[i]:state[i] for i in range(2*self.nmax+1)
            if state[i]!= 0})



def isPinv(state):
    return (state[::-1]==state).all()


def occn(state):
    return sum(state)


class Basis():
    # @profile
    def __init__(self, k, stateList):
        self.k = k

        energy = self.info.energy

        # Order the states in energy
        self.stateList = list(sorted(stateList, key=energy))
        self.size = len(self.stateList)

        self.repr1List = [tuple(self.info.repr2torepr1(state).elements())
                for state in self.stateList]

        # 1 if the state is P invariant, 0 if it's not
        self.parityList = array([int((state == state[::-1]).all())
            for state in self.stateList])
        self.energyList = array([energy(state) for state in self.stateList])
        self.occnList = array([sum(state) for state in self.stateList])

        self.Emax = self.energyList[-1]
        self.Emin = self.energyList[0]

        # Contains also the P-reversed states
        # NOTE: using arrays is much less efficient!
        self.statePos = {}
        for i,state in enumerate(self.stateList):
            self.statePos[tuple(state)] = i
            self.statePos[tuple(state[::-1])] = i
        # self.statePos = {
                # **{tuple(state.tolist()):i for i, state in enumerate(self.stateList)},
                # **{tuple(state[::-1].tolist()):i for i, state in enumerate(self.stateList)}}

        self.stateDlists = {nd: [set(tuple(sorted(x)) for x in combinations(state,nd))
                for state in self.repr1List] for nd in range(5)}


    @classmethod
    def fromScratch(self, Emax, info, occmax=None):
        """ Builds the truncated Hilbert space up to cutoff Emax """

        self.info = info

        self.Emax = Emax
        m = info.m
        L = info.L

        # This can be "None"
        self.occmax = occmax
        if occmax==None:
            # This is used in the basis construction
            self._occmax = int(floor(Emax/m))
        else:
            self._occmax = occmax

        # self.nmax is the actual maximum occupied wavenumber of the states
        self.nmax = int(floor(sqrt((Emax/2.)**2.-m**2.)*L/(2*pi)))

        bases = self.buildBasis(self)
        return {k:self(k,bases[k]) for k in (-1,1)}

    @classmethod
    def fromBasis(self, basis, filterFun):
        """ Extracts a sub-basis with vectors v such that filterFun(v)=True """
        stateList = [v for v in basis.stateList if filterFun(v) == True]
        return self(basis.k, stateList)


    def __repr__(self):
        return str(self.stateList)
    def __getitem__(self,index):
        return self.stateList[index]

    # @profile
    def lookup(self, state):
        """ Looks up the index of a state (array of occupation numbers) """
        # return self.statePos[tuple(state.tolist())]
        x = tuple(state)
        return self.statePos[x]


    def genRMlist(self, RMstate, n):
        """ Recursive function generating all the states starting from RMstate, by adding
        any number of particles with wavenumber n.
        It starts from the seed state with 0 particles and wavenumber 1 """

        omega = self.info.omega

        if n > self.nmax:
            return [RMstate]

        maxN = int(floor(min(
                self._occmax-occn(RMstate),
                (self.nmax-self.info.RMtotalWN(RMstate))/n,
                (self.Emax-self.info.RMenergy(RMstate))/self.info.omega(n))))
        ret = []
        for N in range(maxN+1):
            newstate = scipy.copy(RMstate)
            newstate[n-1] = N
            ret += self.genRMlist(self,newstate,n+1)
        return ret


    def buildBasis(self):
        """ Generates the basis starting from the list of RM states """

        omega = self.info.omega
        m = self.info.m

        RMlist = self.genRMlist(self, scipy.zeros(self.info.nmax,dtype=np.int8), n=1)

        # divides the list of RMstates into a list of lists,
        # so that two states in each list have a fixed total RM wavenumber,
        # and sort each sublist in energy
        sortedRMlist = sorted(RMlist, key=self.info.RMtotalWN)
        dividedRMlist = [sorted(l, key=self.info.RMenergy) for wn,l in
            itertools.groupby(sortedRMlist,key=self.info.RMtotalWN)]
        statelist = {1:[], -1:[]}

        for RMwn, RMsublist in enumerate(dividedRMlist):
            for i, RMstate in enumerate(RMsublist):
                ERM = self.info.RMenergy(RMstate)
                ORM = occn(RMstate)

                # LM part of the state will come from the same sublist.
                # We take the position of LMState to be greater or equal
                # to the position of RMstate
                for LMstate in RMsublist[i:]:
                    # we will just have to reverse it
                    ELM = self.info.RMenergy(LMstate)
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
                        state = scipy.concatenate((LMstate[::-1],array([N0]),RMstate))
                        statelist[(-1)**(N0+OLM+ORM)].append(state)

        return statelist

