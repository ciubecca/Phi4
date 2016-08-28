import scipy
import scipy.sparse.linalg
import scipy.sparse
from scipy import pi
import math
from operator import attrgetter
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
        self.occmax = occmax

        # Maximum wavenumber
        self.nmax = int(math.floor(scipy.sqrt((Emax/2.)**2.-m**2.)*self.L/(2.*pi)))
        if occmax==None:
            self._occmax = int(math.floor(Emax/self.m))
        else:
            self._occmax = occmax

        # Collection of Fock space states, possibly sorted in energy
        stateList = sorted(self.__buildBasis(self), key=attrgetter('energy'))

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


    def __buildRMlist(self):
        """
        sets list of all right -moving states with particles of individual wave number
        <= nmax, total momentum <= Emax/2 and total energy <= Emax
        This function works by first filling in n=1 mode in all possible ways, then n=2 mode
        in all possible ways assuming the occupation of n=1 mode, etc
        """

        kmax = self.nmax*2*pi/self.L

        # maximal occupation number of n=1 mode
        # If there is a right-moving particle, there must be also a left-moving one
        maxN1 = int(math.floor(min(self._occmax-1, kmax/k(1,self.L), self.Emax/omega(1,self.L,self.m))))

        # seed list of RM states, all possible n=1 mode occupation numbers
        RMlist0 = [State([N],1,L=self.L,m=self.m,checkAtRest=False) for N in range(maxN1+1)]

        # go over all other modes
        for n in range(2,self.nmax+1):
            # we will take states out of RMlist0, augment them and add to RMlist1
            RMlist1=[]
            # cycle over all RMstates
            for RMstate in RMlist0:
                p0 = RMstate.momentum
                e0 = RMstate.energy
                # maximal occupation number of mode n given the occupation numbers of
                # all previous modes
                maxNn = int(math.floor(
                    min(self._occmax-RMstate.occ-1, (kmax-p0)/k(n,self.L),
                        (self.Emax-scipy.sqrt(self.m**2+p0**2)-e0)/omega(n,self.L,self.m))))

                for N in range(maxNn+1):
                    longerstate=RMstate.occs[:]
                    #add all possible occupation numbers for mode n
                    longerstate.append(N)
                    RMlist1.append(State(longerstate,len(longerstate),
                        L=self.L,m=self.m, checkAtRest=False))
            #RMlist1 created, copy it back to RMlist0
            RMlist0 = RMlist1

        #save list of RMstates in an internal variable
        self.__RMlist = RMlist0

    def __divideRMlist(self):
        """ divides the list of RMstates into a list of lists, RMdivided,
        so that two states in each list have a fixed total RM wavenumber,
        also each sublist is ordered in energy"""

        self.__nRMmax=max([RMstate.totalWN for RMstate in self.__RMlist])
        #initialize list of lists
        self.__RMdivided = [[] for ntot in range(self.__nRMmax+1)]
        #go over RMstates and append them to corresponding sublists
        for RMstate in self.__RMlist:
            self.__RMdivided[RMstate.totalWN].append(RMstate)

        #now sort each sublist in energy
        for RMsublist in self.__RMdivided:
            RMsublist.sort(key=attrgetter('energy'))

    # finally function which builds the basis
    def __buildBasis(self):
        """
        creates basis of states of total momentum zero, energy <= Emax
        and occupation number <= occmax
        """
        self.__buildRMlist(self)
        self.__divideRMlist(self)

        statelist = []

        for nRM,RMsublist in enumerate(self.__RMdivided):
            for i, RMstate in enumerate(RMsublist):
                ERM = RMstate.energy
                ORM = RMstate.occ

                # LM part of the state will come from the same sublist.
                # We take the position of LMState to be greater or equal
                #to the position of RMstate
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

                        state = State(LMstate.occs[::-1]+[N0]+RMstate.occs, self.nmax, L=self.L,m=self.m,checkAtRest=True)

                        # XXX just to be sure
                        if state.energy <= self.Emax and state.occ <= self._occmax:
                            # XXX Can this be optimized?
                            if self.k == state.kparity():
                                statelist.append(state)
        return statelist
