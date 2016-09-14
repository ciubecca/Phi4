import scipy
import scipy.sparse.linalg
import scipy.sparse
from scipy import pi
import math
from operator import attrgetter
import itertools
import gc

# "tolerance" parameter to avoid possible numerical
# issues when confronting energies with Emax
tol = 10**(-5)

"""
P denotes spatial parity, while K field parity.
For now only the P-even sector is implemented
"""

def omega(n,L,m):
    """ computes one particle energy from wavenumber"""
    return math.sqrt(m**2.+((2.*pi/L)*n)**2.)
def k(n,L):
    """ computes momentum from wavenumber"""
    return (2.*pi/L)*n

class State():
    def __init__(self, occs, nmax, L=None, m=None, fast=False, checkAtRest=True):
        """ occs: occupation number list
            nmax: wave number of the last element in occs """
        self.occs = occs
        self.nmax = nmax
        self.fast = fast
        self.size = len(self.occs)

        if(fast == False):
            self.setLM(L,m, checkAtRest)

    def setLM(self, L, m, checkAtRest=True):
        self.fast = False

        wavenum = scipy.vectorize(lambda i: i-self.size+self.nmax+1)(range(self.size))
        energies = scipy.vectorize(lambda k : omega(k,L,m))(wavenum)

        self.totalWN = (wavenum*self.occs).sum()

        if checkAtRest and self.totalWN != 0:
            raise ValueError("State not at rest")

        if self.size == 2*self.nmax+1 and self.occs[::-1] == self.occs:
            self.__parityEigenstate = True
        else:
            self.__parityEigenstate = False

        self.L=L
        self.m=m
        self.energy = sum(energies*self.occs)
        # Total occupation number
        self.occ = sum(self.occs)
        self.momentum = (2.*pi/self.L)*self.totalWN

    def isParityEigenstate(self):
        """ Returns True if the Fock space state is a P-parity eigenstate """
        return self.__parityEigenstate

    def kparity(self):
        """ Returns the K-parity quantum number """
        return (-1)**sum(self.occs)

    def __repr__(self):
        return str(self.occs)
    def __eq__(self, other):
       # check also if the P-reversed is the same!
       # return (self.occs == other.occs).all() or (self.occs == other.occs[::-1]).all()
       return (self.occs == other.occs) or (self.occs == other.occs[::-1])

    def __hash__(self):
        return hash(tuple(self.occs))

    def __setitem__(self, wn, n):
        """ Sets the occupation number corresponding to a wave number """
        if self.fast == False:
            self.energy += (n-self[wn])*omega(wn,self.L,self.m)
            self.totalWN += (n-self[wn])*wn
            self.momentum = (2.*pi/self.L)*self.totalWN

        self.occs[wn+self.size-self.nmax-1] = n

    def __getitem__(self, wn):
        """ Returns the occupation number corresponding to a wave number """
        return self.occs[wn+self.size-self.nmax-1]

    def wnList(self):
        return range(-self.nmax,self.nmax+1)

    def parityReversed(self):
        """ Reverse under P parity """
        if not self.size == 2*self.nmax+1:
            raise ValueError("attempt to reverse asymmetric occupation list")
        return State(self.occs[::-1],self.nmax,L=self.L,m=self.m)


class Basis():
    def __init__(self, L, m, k, stateList, Emax=None):
        self.L = L
        self.Emax = Emax
        self.m = m
        self.k = k

        self.stateList = stateList
        self.size = len(self.stateList)

        self.statePos = {state:i for i, state in enumerate(self.stateList) }

    @classmethod
    def fromScratch(self, L, Emax, m, k, occmax=None):
        """ Builds the truncated Hilbert space up to cutoff Emax """
        self.L = L
        self.Emax = Emax
        self.m = m
        self.k = k
        # This can be "None"
        self.occmax = occmax

        if occmax==None:
            # This is used in the basis construction
            self._occmax = int(math.floor(Emax/self.m))
        else:
            self._occmax = occmax

        self.kmax = scipy.sqrt((Emax/2.)**2.-m**2.)
        self.nmax = int(math.floor(self.kmax*L/(2*pi)))

        # Collection of Fock space states, possibly sorted in energy
        stateList = sorted(self.buildBasis(self), key=attrgetter('energy'))

        return self(L, m, k, stateList, Emax)

    @classmethod
    def fromBasis(self, basis, filterFun):
        """ Extracts a sub-basis with vectors v such that filterFun(v)=True """
        stateList = [v for v in basis.stateList if filterFun(v) == True]
        return self(basis.L, basis.m, basis.k, stateList)


    def __len__(self):
        return len(self.stateList)
    def __repr__(self):
        return str(self.stateList)
    def __getitem__(self,index):
        return self.stateList[index]

    def lookup(self, state):
        """ Looks up the index of a state. """
        try:
            return self.statePos[state]
        except KeyError:
            raise LookupError()


    def genRMlist(self, RMstate, n):
        """ Recursive function generating all the states starting from RMstate, by adding
        any number of particles with wavenumber n.
        It starts from the seed state with 0 particles and wavenumber 1 """

        if n > self.nmax:
            return [RMstate]

        maxN = int(math.floor(min(
                self._occmax-RMstate.occ,
                (self.kmax-RMstate.momentum)/k(n,self.L),
                (self.Emax-RMstate.energy)/omega(n,self.L,self.m))))

        ret = []
        for N in range(maxN+1):
            newoccs = RMstate.occs[:]
            newoccs[n-1] = N
            newstate = State(newoccs, self.nmax, L=self.L, m=self.m, checkAtRest=False)
            ret += self.genRMlist(self,newstate,n+1)
        return ret


    def buildBasis(self):
        """ Generates the basis starting from the list of RM states """

        RMlist = self.genRMlist(self,
                RMstate=State([0]*self.nmax,self.nmax,self.L,self.m), n=1)

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
                ORM = RMstate.occ

                # LM part of the state will come from the same sublist.
                # We take the position of LMState to be greater or equal
                # to the position of RMstate
                for LMstate in RMsublist[i:]:
                    # we will just have to reverse it
                    ELM = LMstate.energy
                    OLM = LMstate.occ
                    deltaE = self.Emax - ERM - ELM
                    deltaOcc = self._occmax - ORM - OLM

                    # if this happens, we can break since subsequent LMstates
                    # have even higher energy (RMsublist is ordered in energy)
                    if deltaE<0:
                        break
                    if deltaOcc<0:
                        continue

                    maxN0 = min(deltaOcc, math.floor(deltaE/self.m))

                    # possible values for the occupation value at rest
                    for N0 in range(maxN0+1):
                        # Only states with correct parity
                        if (-1)**(N0+OLM+ORM) == self.k:
                            statelist.append(State(LMstate.occs[::-1]+[N0]+RMstate.occs,
                                        self.nmax, L=self.L, m=self.m, checkAtRest=True))

        return statelist
