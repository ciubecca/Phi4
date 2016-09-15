import scipy
from math import factorial, floor, sqrt
from statefuncs import Basis, omega, State, k
from collections import Counter
from itertools import combinations
from operator import attrgetter
from matrix import Matrix
from scipy import exp, pi, array
from sortedcontainers import SortedList


def energy(wnlist):
    return sum(omega(k(n)) for n in wnlist)


def bose(x):
    """ computes the Bose factor of a product of oscillators  """
    return factorial(len(x))/scipy.prod(list(map(factorial,Counter(x).values())))

class Operator():
    """
    Collection of oscillators with fixed number of creation and annihilation operators
    """

    def __init__(self, oscillators, nd, nc):
        """
        oscillators: dict where each key is a SORTED tuple of wavenumbers of
        annihilation operators, and each value is a list of tuples of wavenumbers
        of creation operators
        nd: number of annihilation operators
        nc: number of creation operators
        """

        self.nd = nd
        self.nc = nc

        # Sort the lists of creation operator in energy for performance reasons
        self.oscillators = {dlist: list(sorted(clists,key=energy))
            for dlist,clists in oscillators.items()}

        # Define auxiliary structures
        self.oscEnergies = {dlist: list(map(lambda x: -energy(dlist)+energy(x), clists))
            for dlist,clists in self.oscillators.items()}


        # Multiply by Bose factors, the factor coming normal ordering
        # and "phase space" factor
        self.oscFactors = {dlist:
                list(map(lambda x: bose(x)*bose(dlist)*scipy.special.binom(nc+nd,nc)\
            *scipy.prod([1/sqrt(2*omega(k(n))*L) for n in x+dlist]), clists))
            for dlist,clists in self.oscillators.items()}


        #XXX Define oscillators in the "1" representation


    @profile
    def computeMatrixElements(self, state, lookupbasis):
        # List of columns indices of generated basis elements
        col = []
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

            oscEnergies = self.oscEnergies[dlist]
            oscFactors = self.oscFactors[dlist]

            for i,clist in enumerate(self.oscillators[dlist]):
                # XXX Can be optimized with SortedList
                if e + oscEnergies[i] > Emax:
                    break

                # print("dlist", dlist)
                # print("clist", clist)
                # print("e", e)
                # print("oscE", oscEnergies[i])


                coeff = oscFactors[i]
                # print(coeff)

                # XXX Can this be optimized using a table and a difference of vectors?
                newstate = occs.tolist()

                for d in dlist:
                    coeff *= sqrt(newstate[d+nmax])
                    newstate[d+nmax] -= 1
                for c in clist:
                    coeff *= sqrt(newstate[c+nmax]+1)
                    newstate[c+nmax] += 1

                j = lookup(newstate)
                # j = lookupdict[tuple(newstate)]

                # parity factors XXX check
                if p:
                    coeff *= 1/sqrt(2)
                if lookupbasis[j].isPeigenstate:
                    coeff *= sqrt(2)


                # print(coeff)
                col.append(j)
                data.append(coeff)

        return col, data

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

    V40 = Operator(V40, 0, 4)


    V31 = {}
    for k1 in range(-nmax,nmax+1):
        dlist = (k1,)
        V31[dlist] = []

        for k2 in range(-nmax,nmax+1):
            for k3 in range(max(-nmax+k1-k2,k2),
                           min(int(floor((k1-k2)/2)),nmax)+1):

                k4 = k1-k2-k3
                clist = (k2,k3,k4)

                if -Emax <= energy(clist)-energy(dlist) <= Emax:
                    V31[dlist].append(clist)

    V31 = Operator(V31, 1, 3)


    V22 = {}
    for k1 in range(-nmax,nmax+1):
        for k2 in range(k1,nmax+1):
            dlist = (k1,k2)
            V22[dlist] = []

            for k3 in range(max(-nmax+k1+k2,-nmax),
                    min(int(floor((k1+k2)/2)),nmax)+1):

                k4 = k1+k2-k3
                clist = (k3,k4)

                if -Emax <= energy(clist)-energy(dlist) <= Emax and\
                    sorted([abs(k3),abs(k4)])<=sorted([abs(k1),abs(k2)]):
                    # only consider lexicographically ordered part of V22
                    # but also including diagonal part which will be separated below

                    V22[dlist].append(clist)

    V22 = Operator(V22, 2, 2)

    return V40, V31, V22
