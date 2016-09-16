import scipy
from math import factorial, floor, sqrt
from statefuncs import Basis, omega, State, k
from collections import Counter
from itertools import combinations
from operator import attrgetter
from matrix import Matrix
from scipy import exp, pi, array
from sortedcontainers import SortedList
import bisect

def energy(wnlist):
    return sum(omega(k(n)) for n in wnlist)


def bose(x):
    """ computes the Bose factor of a product of oscillators  """
    return factorial(len(x))/scipy.prod(list(map(factorial,Counter(x).values())))


def repr1torepr2(state):
    """ Transform state from repr1 to repr2 """
    return array([state.get(n,0) for n in range(-nmax,nmax+1)])


def subtractRepr3toRepr1(osc1, osc2):
    ret = Counter(osc1)
    ret.subtract(Counter(osc2))
    return Counter({n:Zn for n,Zn in ret.items() if Zn!=0})

class Operator():
    """
    Collection of oscillators with fixed number of creation and annihilation operators
    """

    def __init__(self, oscillators, nd, nc, basis):
        """
        oscillators: dict where each key is a SORTED tuple of wavenumbers of
        annihilation operators, and each value is a list of tuples of wavenumbers
        of creation operators
        nd: number of annihilation operators
        nc: number of creation operators
        """

        self.nd = nd
        self.nc = nc

        # TODO we could define just one dict of indices instead of many dicts

        # Sort the lists of creation operator in energy for performance reasons
        self.oscillators = {dlist: list(sorted(clists,key=energy))
            for dlist,clists in oscillators.items()}


        # Store the values of the energies of the oscillators
        self.oscEnergies = {dlist:
                list(map(lambda clist: -energy(dlist)+energy(clist), clists))
            for dlist,clists in self.oscillators.items()}


        # Multiply by Bose factors, the factor coming normal ordering
        # and "phase space" factor
        self.oscFactors = {dlist:
                list(map(lambda clist: bose(clist)*bose(dlist)\
                    *scipy.special.binom(nc+nd,nc)\
                    *scipy.prod([1/sqrt(2*omega(k(n))*L) for n in clist+dlist]), clists))
            for dlist,clists in self.oscillators.items()}

#         # Store the "difference" states of the oscillators in repr "1"
        # self.doscRepr1 = {dlist: Counter(dlist) for dlist in self.oscillators.keys()}

        # # Store the "difference" states of the oscillators in repr "1"
        # self.coscRepr1 = {dlist:
                # list(map(Counter, clists))
            # for dlist,clists in self.oscillators.items()}

        self.doscRepr2 = {dlist: repr1torepr2(Counter(dlist))
                for dlist in self.oscillators.keys()}

        self.coscRepr2 = {dlist: array(list(map(lambda clist:
            repr1torepr2(Counter(clist)), clists)))
                for dlist,clists in self.oscillators.items()}



        # Store the "difference" states of the oscillators in repr "1"
        self.oscRepr1 = {dlist:
                list(map(lambda clist: subtractRepr3toRepr1(clist,dlist), clists))
           for dlist,clists in self.oscillators.items()}


        # Store the "difference" states of the oscillators in repr "2"
        self.oscRepr2 = {dlist:
                array(list(map(repr1torepr2, states))) for dlist,states in self.oscRepr1.items()}



        # Precompute normalization factors
        maxocc = 30
        # c, d, n
        self.normFactors = scipy.zeros(shape=(nc+1,nd+1,maxocc+1))
        for c in range(nc+1):
            for d in range(nd+1):
                for n in range(maxocc+1):
                    if d <= n:
                        self.normFactors[c,d,n] = \
                            sqrt(factorial(n)*factorial(n-d+c)/factorial(n-d)**2)
                    else:
                        self.normFactors[c,d,n] = scipy.nan

        # Precompute parity factors
        # XXX check
        self.parityFactors = scipy.array([[1,sqrt(2)],[1/sqrt(2),1]])

        self.P = array([int(state.isPeigenstate) for state in basis.stateList])

    # @profile
    def computeMatrixElements(self, state, lookupbasis):
        # List of columns indices of generated basis elements
        col = []
        col2 = []
        data2 = []
        # List of partial matrix elements
        data = []

        # I define these local variables for performance reasons
        p = state.isPeigenstate
        e = state.energy
        Emax = lookupbasis.Emax
        occs = state.occs
        lookup = lookupbasis.lookup

        # print("state", state)
        # print("state wnlist", state.wnlist)

        # Select all the possible groups of particles in the state than can be annihilated
        # Call to set() is needed to eliminate duplicates
        for dlist in set(combinations(state.wnlist, self.nd)):


            # print("dlist", dlist)

            imax = bisect.bisect_left(self.oscEnergies[dlist], Emax-e+(10**-10))
            if imax==0:
                continue

            oscFactors = self.oscFactors[dlist][:imax]
            oscRepr2 = self.oscRepr2[dlist][:imax]
            # oscillators = self.oscillators[dlist][:imax]
            # coscRepr1 = self.coscRepr1[dlist][:imax]
            # doscRepr1 = self.doscRepr1[dlist]
            coscRepr2 = self.coscRepr2[dlist][:imax]
            doscRepr2 = self.doscRepr2[dlist]


            newstateList = oscRepr2 + occs
            colpart = array([lookup(state) for state in newstateList])
            # colpart = lookup(newstateList)
            # print("nmax", nmax)
            # print("n osc", len(oscRepr2))
            # print("data shape", fv(occs, coscRepr2, doscRepr2).shape)
            # print(scipy.prod(fv(occs, coscRepr2, doscRepr2)))
            normFactors = self.normFactors[coscRepr2,doscRepr2,newstateList]
            datapart = scipy.prod(normFactors, axis=1).tolist()
            # data2 += scipy.prod(fv(occs, coscRepr2, doscRepr2), axis=1).tolist()
            # print(data2.shape)

            datapart *= self.parityFactors[int(p),self.P[colpart]]

            col2 += colpart.tolist()
            data2 += datapart.tolist()

#             for i,clist in enumerate(oscillators):
                # # XXX Can be optimized with SortedList
                # # if e + oscEnergies[i] > Emax:
                    # # break

                # # print("dlist", dlist)
                # # print("clist", clist)
                # # print("e", e)
                # # print("oscE", oscEnergies[i])


                # coeff = oscFactors[i]
                # # print(coeff)

                # newstate2 = occs-oscRepr2[i]


                # # XXX Can this be optimized using a table and a difference of vectors?
                # newstate = occs.tolist()
                # for d in dlist:
                    # coeff *= sqrt(newstate[d+nmax])
                    # newstate[d+nmax] -= 1
                # for c in clist:
                    # coeff *= sqrt(newstate[c+nmax]+1)
                    # newstate[c+nmax] += 1

                # j = lookup(newstate)
                # # j = lookupdict[tuple(newstate)]

                # # parity factors XXX check
                # if p:
                    # coeff *= 1/sqrt(2)
                # if lookupbasis[j].isPeigenstate:
                    # coeff *= sqrt(2)


                # # print(coeff)
                # col.append(j)
                # data.append(coeff)

        # print(len(data2)==len(data))

        return col2, data2

    # TODO Generate high energy Hilbert space Hh from low energy Hilbert space Hl
    # as Hh = V*Hl
    def genBasis(self):
        return


def Phi4Operators(lookupbasis, LL, mm, nmaxx):

    global L, m, nmax
    L = LL
    m = mm
    nmax = nmaxx

    Emax = lookupbasis.Emax

    dlist = ()
    V40 = {dlist: []}
    for k1 in range(-nmax,nmax+1):
        for k2 in range(k1,nmax+1):
            # NOTE the boundaries for k3 ensure that k3<=k4<=nmax
            for k3 in range(max(-nmax-k1-k2,k2),
                    min(int(floor((-k1-k2)/2)),nmax)+1):

                k4 = -k1-k2-k3
                clist = (k1,k2,k3,k4)

                if energy(clist) <= Emax:
                    V40[dlist].append(clist)

    V40 = Operator(V40, 0, 4, lookupbasis)


    V31 = {}
    for k1 in range(-nmax,nmax+1):
        dlist = (k1,)
        V31[dlist] = []

        for k2 in range(-nmax,nmax+1):
            for k3 in range(max(-nmax+k1-k2,k2),
                           min(int(floor((k1-k2)/2)),nmax)+1):

                k4 = k1-k2-k3
                clist = (k2,k3,k4)

                # NOTE The check on dlist is useless here but it'ŝ needed
                # if we generalize the code to other operators
                if energy(clist) <= Emax and energy(dlist) <= Emax:
                    V31[dlist].append(clist)

    V31 = Operator(V31, 1, 3, lookupbasis)


    V22 = {}
    for k1 in range(-nmax,nmax+1):
        for k2 in range(k1,nmax+1):
            dlist = (k1,k2)
            V22[dlist] = []

            for k3 in range(max(-nmax+k1+k2,-nmax),
                    min(int(floor((k1+k2)/2)),nmax)+1):

                k4 = k1+k2-k3
                clist = (k3,k4)

                if energy(dlist) <= Emax and energy(clist) <= Emax and\
                    sorted([abs(k3),abs(k4)])<=sorted([abs(k1),abs(k2)]):
                    # only consider lexicographically ordered part of V22
                    # but also including diagonal part which will be separated below

                    V22[dlist].append(clist)

    V22 = Operator(V22, 2, 2, lookupbasis)

    return V40, V31, V22
