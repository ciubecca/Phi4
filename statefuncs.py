import scipy
from scipy import array, pi, sqrt
import math
from operator import attrgetter
import itertools
import numpy as np

# Global variables
L = None
m = None

def omega(k):
    """ computes one particle energy from momentum"""
    return sqrt(m**2.+k**2.)
def k(n):
    """ computes momentum from wavenumber"""
    return (2*pi/L)*n


class State():
    def __init__(self, occs, nmax, atRest=True):
        """ occs: occupation number list
            nmax: wave number of the last element in occs """
        self.occs = occs

        self.nmax = nmax
        wavenums = array(range(self.occs.size))-self.occs.size+self.nmax+1
        self.energy = (self.occs*omega(k(wavenums))).sum()
        self.totalWN = (wavenums*self.occs).sum()
        self.occn = self.occs.sum()

        if atRest:
            if self.totalWN != 0:
                raise ValueError("State not at rest")

            self.kparity = (-1)**self.occn
            self.isPeigenstate = (self.occs == self.occs[::-1]).all()

            # Alternative representation of the state: ordered list of
            # occupied wavenumbers (can contain duplicates)
            self.wnlist = []
            for n in wavenums:
                self.wnlist += [n]*self[n]


    def __repr__(self):
        return str(self.occs)


    def __getitem__(self, wn):
        """ Returns the occupation number corresponding to a wave number """
        return self.occs[wn+self.occs.size-self.nmax-1]

    def Preversed(self):
        return State(self.occs[::-1],self.nmax,atRest=False)


class Basis():
    def __init__(self, k, stateList, Emax=None):
        self.Emax = Emax
        self.k = k

        self.stateList = stateList
        self.size = len(self.stateList)

        # Contains also the P-reversed states
        # NOTE: For some reason using arrays is much less efficient!
        self.statePos = {
                **{tuple(state.occs.tolist()):i for i, state in enumerate(self.stateList)},
                **{tuple(state.Preversed().occs.tolist()):i for i, state in
                        enumerate(self.stateList) if not state.isPeigenstate}
                }

    @classmethod
    def fromScratch(self, LL, mm, Emax, k, occmax=None):
        """ Builds the truncated Hilbert space up to cutoff Emax """
        global L
        L = LL
        global m
        m = mm

        self.Emax = Emax
        self.k = k
        # This can be "None"
        self.occmax = occmax

        if occmax==None:
            # This is used in the basis construction
            self._occmax = int(math.floor(Emax/m))
        else:
            self._occmax = occmax

        kmax = scipy.sqrt((Emax/2.)**2.-m**2.)
        self.nmax = int(math.floor(kmax*L/(2*pi)))

        # Collection of Fock space states, possibly sorted in energy
        stateList = sorted(self.buildBasis(self), key=attrgetter('energy'))

        return self(k, stateList, Emax)

    @classmethod
    def fromBasis(self, basis, filterFun):
        """ Extracts a sub-basis with vectors v such that filterFun(v)=True """
        stateList = [v for v in basis.stateList if filterFun(v) == True]
        return self(basis.k, stateList)


    def __len__(self):
        return len(self.stateList)
    def __repr__(self):
        return str(self.stateList)
    def __getitem__(self,index):
        return self.stateList[index]

    def lookup(self, state):
        """ Looks up the index of a state (list of occupation numbers) """
        return self.statePos[tuple(state)]


    def genRMlist(self, RMstate, n):
        """ Recursive function generating all the states starting from RMstate, by adding
        any number of particles with wavenumber n.
        It starts from the seed state with 0 particles and wavenumber 1 """


        if n > self.nmax:
            return [RMstate]

        maxN = int(math.floor(min(
                self._occmax-RMstate.occn,
                (self.nmax-RMstate.totalWN)/n,
                (self.Emax-RMstate.energy)/omega(k(n)))))
        ret = []
        for N in range(maxN+1):
            newoccs = scipy.copy(RMstate.occs)
            newoccs[n-1] = N
            newstate = State(newoccs, self.nmax, atRest=False)
            ret += self.genRMlist(self,newstate,n+1)
        return ret


    def buildBasis(self):
        """ Generates the basis starting from the list of RM states """

        RMlist = self.genRMlist(self,
                RMstate=State(scipy.zeros(self.nmax,dtype=np.int8),self.nmax,atRest=False),
                n=1)

        # divides the list of RMstates into a list of lists,
        # so that two states in each list have a fixed total RM wavenumber,
        # and sort each sublist in energy
        sortedRMlist = sorted(RMlist, key=attrgetter('totalWN'))
        dividedRMlist = [sorted(l, key=attrgetter('energy')) for wn,l in
            itertools.groupby(sortedRMlist,key=attrgetter('totalWN'))]
        statelist = []

        for RMwn, RMsublist in enumerate(dividedRMlist):
            for i, RMstate in enumerate(RMsublist):
                ERM = RMstate.energy
                ORM = RMstate.occn

                # LM part of the state will come from the same sublist.
                # We take the position of LMState to be greater or equal
                # to the position of RMstate
                for LMstate in RMsublist[i:]:
                    # we will just have to reverse it
                    ELM = LMstate.energy
                    OLM = LMstate.occn
                    deltaE = self.Emax - ERM - ELM
                    deltaOcc = self._occmax - ORM - OLM

                    # if this happens, we can break since subsequent LMstates
                    # have even higher energy (RMsublist is ordered in energy)
                    if deltaE<0:
                        break
                    if deltaOcc<0:
                        continue

                    maxN0 = int(math.floor(min(deltaOcc, deltaE/m)))

                    # possible values for the occupation value at rest
                    for N0 in range(maxN0+1):
                        # Only states with correct parity
                        if (-1)**(N0+OLM+ORM) == self.k:
                            state = scipy.concatenate((LMstate.occs[::-1],array([N0]),
                                    RMstate.occs))
                            statelist.append(State(state, self.nmax))

        return statelist
