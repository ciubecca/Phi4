import scipy
from math import factorial, floor, sqrt
from statefuncs import Basis
from collections import Counter
import itertools
from itertools import combinations, islice
from scipy import exp, pi, array
import bisect

tol = 10**(-10)

# @profile
def gendlists(state, nd):
    x = list(itertools.chain.from_iterable([[n]*Zn for n,Zn in state]))
    if nd==4:
        ret = set(tuple(y) for y in combinations(x,nd))
        return set(dlist for dlist in ret if sum(dlist)==0)
    else:
        ret = set(tuple(y) for y in combinations(x,nd))
        return ret

def bose(x):
    """ computes the Bose factor of a product of oscillators  """
    return factorial(len(x))/scipy.prod(list(map(factorial,Counter(x).values())))


# XXX Check
parityFactors = [[1, sqrt(2)],[1/sqrt(2),1]]


class Operator():
    """
    Collection of oscillators with fixed number of creation and annihilation operators
    """

    def __init__(self, oscillators, nd, nc, helper):
        """
        oscillators: list of tuples. The first element of the tuple is a tuple of
        wavenumbers of annihilation operators, and the second element a list of
        tuples of wavenumbers of creation operators
        basis: basis on which the operator will act
        nd: number of annihilation operators
        nc: number of creation operators
        """

        self.nd = nd
        self.nc = nc
        omega = helper.omega
        L = helper.L
        m = helper.m
        self.helper = helper

        self.dlistPos = {}
        self.oscList = []
        self.oscEnergies = []
        self.oscFactors = []

        def f(clist, dlist):
            ret = []
            wnlist = set(clist+dlist)
            cosc = Counter(clist)
            dosc = Counter(dlist)
            for n in wnlist:
                ret.append((n, cosc.get(n,0), dosc.get(n,0)))
            return ret


        for i, (dlist,clists) in enumerate(oscillators):
            clists = list(sorted(clists,key=helper.oscEnergy))

            self.dlistPos[dlist] = i

            self.oscList.append([f(clist,dlist) for clist in clists])

            self.oscEnergies.append([helper.oscEnergy(clist)-helper.oscEnergy(dlist)
                for clist in clists])

            self.oscFactors.append([bose(clist)*bose(dlist)*scipy.special.binom(nc+nd,nc)\
                    *scipy.prod([1/sqrt(2*omega(n)*L) for n in clist+dlist])
                    for clist in clists])


    # @profile
    def computeMatrixElements(self, basis, i, lookupbasis, statePos=None):
        # List of columns indices of generated basis elements
        col = []
        # List of partial matrix elements
        data = []


        # I define these local variables outside the loops for performance reasons
        e = basis.energyList[i]
        p = basis.parityList[i]
        state = basis.stateList[i]
        # statevec = self.helper.torepr2(state)
        if basis.helper.nmax > lookupbasis.helper.nmax:
            helper = basis.helper
        else:
            helper = lookupbasis.helper

        statevec = helper.torepr2(state)
        parityList = lookupbasis.parityList
        Emin = lookupbasis.Emin
        Emax = lookupbasis.Emax
        nmax = helper.nmax


        if statePos == None:
            statePos = lookupbasis.statePos


        normFactors = self.helper.normFactors

        # cycle over all the sets of momenta that can be annihilated
        for dlist in gendlists(state, self.nd):

            k = self.dlistPos[dlist]

            imin = bisect.bisect_left(self.oscEnergies[k], Emin-e-tol)
            imax = bisect.bisect_left(self.oscEnergies[k], Emax-e+tol)

            if imax <= imin:
                continue

            datapart = (self.oscFactors[k][imin:imax])[:]

            # This cycle could be moved to C
            for i, osc in enumerate(self.oscList[k][imin:imax]):
                newstatevec = statevec[:]

                for n,Zc,Zd in osc:
                    newstatevec[n+nmax] += Zc-Zd
                    datapart[i] *= normFactors[Zc, Zd, statevec[n+nmax]]
                try:
                    j = statePos[tuple(newstatevec)]
                except KeyError as err:
                    try:
                        newstate = helper.torepr1(newstatevec)
                        print(newstate)
                        print(helper.energy(newstate))
                        raise err
                    except IndexError as err:
                        newstate = helper.torepr1(newstatevec)
                        print("oscillator", osc)
                        # print("osc energy", basis.helper.oscEnergy(osc))
                        print("state", state)
                        print("newstate", newstate)
                        print("state energy", helper.energy(state))
                        print("newstate energy", helper.energy(newstate))
                        raise err

                datapart[i] *= parityFactors[p][parityList[j]]
                col.append(j)

            data += datapart

        return col, data


    # Generate high energy Hilbert space Hh from low energy Hilbert space Hl
    # as Hh = V*Hl
    def genBasis(self, basis, Emin, Emax):

        nmax = self.helper.nmax
        stateset = set()

        for i, state in enumerate(basis):

            statevec = self.helper.torepr2(state)
            e = basis.energyList[i]

            for dlist in gendlists(state, self.nd):
                k = self.dlistPos[dlist]

                # XXX Check
                imin = bisect.bisect_left(self.oscEnergies[k], Emin-e+tol)
                imax = bisect.bisect_left(self.oscEnergies[k], Emax-e+tol)
                if imax <= imin:
                    continue

                for i, osc in enumerate(self.oscList[k][imin:imax]):
                    newstatevec = statevec[:]
                    for n,Zc,Zd in osc:
                        newstatevec[n+nmax] += Zc-Zd
                    t1 = tuple(newstatevec)
                    t2 = tuple(newstatevec[::-1])
                    if (t1 not in stateset) and (t2 not in stateset):
                        stateset.add(t1)

        return stateset



def V4Operators(helper, basis):

    nmax = max(basis[1].nmax, basis[-1].nmax)
    Emax = max(basis[1].Emax, basis[-1].Emax)

    dlist = ()
    V40 = [(dlist, [])]
    for k1 in range(-nmax,nmax+1):
        for k2 in range(k1,nmax+1):
            # NOTE the boundaries for k3 ensure that k3<=k4<=nmax
            for k3 in range(max(-nmax-k1-k2,k2),
                    min(int(floor((-k1-k2)/2)),nmax)+1):

                k4 = -k1-k2-k3
                clist = (k1,k2,k3,k4)

                if helper.oscEnergy(clist) <= Emax+tol:
                    V40[-1][1].append(clist)

    V40 = Operator(V40, 0, 4, helper)


    V31 = []
    for k1 in range(-nmax,nmax+1):
        dlist = (k1,)
        V31.append((dlist,[]))

        for k2 in range(-nmax,nmax+1):
            for k3 in range(max(-nmax+k1-k2,k2),
                           min(int(floor((k1-k2)/2)),nmax)+1):

                k4 = k1-k2-k3
                clist = (k2,k3,k4)
                # NOTE The check on dlist is useless here but it'ŝ needed
                # if we generalize the code to other operators
                if helper.oscEnergy(clist) <= Emax+tol\
                    and helper.oscEnergy(dlist) <= Emax+tol:
                    V31[-1][1].append(clist)

    V31 = Operator(V31, 1, 3, helper)


    V22 = []
    for k1 in range(-nmax,nmax+1):
        for k2 in range(k1,nmax+1):
            dlist = (k1,k2)
            V22.append((dlist,[]))

            for k3 in range(max(-nmax+k1+k2,-nmax),
                    min(int(floor((k1+k2)/2)),nmax)+1):

                k4 = k1+k2-k3
                clist = (k3,k4)

                # NOTE The check on dlist is useless here but it'ŝ needed
                # if we generalize the code to other operators
                if helper.oscEnergy(dlist) <= Emax+tol and\
                    helper.oscEnergy(clist) <= Emax+tol and\
                    sorted([abs(k3),abs(k4)])<=sorted([abs(k1),abs(k2)]):
                    # only consider lexicographically ordered part of V22
                    # but also including diagonal part which will be separated below

                    V22[-1][1].append(clist)


    V22 = Operator(V22, 2, 2, helper)

    return V40, V31, V22



def V2Operators(helper, basis):
    """
    basis: set of states on which the operators act
    """
    nmax = max(basis[1].nmax, basis[-1].nmax)
    Emax = max(basis[1].Emax, basis[-1].Emax)

    dlist = ()
    V20 = [(dlist, [])]
    for k1 in range(-nmax,1):
        k2 = -k1
        clist = (k1,k2)

        # NOTE not needed
        if helper.oscEnergy(clist) <= Emax+tol:
            V20[-1][1].append(clist)

    V20 = Operator(V20, 0, 2, helper)


    V11 = []
    for k1 in range(-nmax,nmax+1):
        dlist = (k1,)
        V11.append((dlist,[]))

        k2 = k1
        clist = (k2,)
        # NOTE Not needed
        if helper.oscEnergy(clist) <= Emax+tol and helper.oscEnergy(dlist) <= Emax+tol:
            V11[-1][1].append(clist)

    V11 = Operator(V11, 1, 1, helper)

    return V20, V11


def V4OperatorsLH(helper, basis, Emin, Emax):
    """
    basis: set of states on which the operators act
    Emin: minimum energy of the states corresponding to columns indices of the matrix
    Emax: maximal energy of the states corresponding to columns indices of the matrix
    """
    nmax = helper.nmax

    dlist = ()
    V40 = [(dlist, [])]
    for k1 in range(-nmax,nmax+1):
        for k2 in range(k1,nmax+1):
            # NOTE the boundaries for k3 ensure that k3<=k4<=nmax
            for k3 in range(max(-nmax-k1-k2,k2),
                    min(int(floor((-k1-k2)/2)),nmax)+1):

                k4 = -k1-k2-k3
                clist = (k1,k2,k3,k4)

                if helper.oscEnergy(clist) <= Emax+tol:
                    V40[-1][1].append(clist)

    V40 = Operator(V40, 0, 4, helper)


    V31 = []

    dlists = set()
    for state in basis.stateList:
        dlists.update(gendlists(state, 1))

    for dlist in dlists:
        k1 = dlist[0]
        V31.append((dlist,[]))

        for k2 in range(-nmax,nmax+1):
            for k3 in range(max(-nmax+k1-k2,k2),
                           min(int(floor((k1-k2)/2)),nmax)+1):

                k4 = k1-k2-k3
                clist = (k2,k3,k4)

                if helper.oscEnergy(clist) <= Emax+tol\
                    and helper.oscEnergy(clist)+tol > helper.oscEnergy(dlist):
                    # We always want to increase the energy
                    V31[-1][1].append(clist)

    V31 = Operator(V31, 1, 3, helper)

    V22 = []

    dlists = set()
    for state in basis.stateList:
        dlists.update(gendlists(state, 2))

    for dlist in dlists:
        (k1,k2) = dlist
        V22.append((dlist,[]))

        for k3 in range(max(-nmax+k1+k2,-nmax),
                min(int(floor((k1+k2)/2)),nmax)+1):

            k4 = k1+k2-k3
            clist = (k3,k4)

            # We always want to increase the energy
            if helper.oscEnergy(clist) >= helper.oscEnergy(dlist) and\
                    helper.oscEnergy(clist) <= Emax+tol:

                V22[-1][1].append(clist)


    V22 = Operator(V22, 2, 2, helper)


    # NOTE V13 and V04 cannot increase the energy

    return V40, V31, V22


# Generate all the oscillators between the "selected" high-energy basis and the
# low-energy basis
def V4OperatorsHL(helper, basis, Emax):
    """
    basis: set of states on which the operators act
    Emax: maximal energy of the states corresponding to columns indices of the matrix
    """

    nmax = helper.nmax

    V04 = []

    dlists = set()
    for state in basis.stateList:
    # This has to be implemented differently
        dlists.update(gendlists(state, 4))

    for dlist in dlists:
        V04.append((dlist,[()]))

    V04 = Operator(V04, 4, 0, helper)


    V13 = []

    dlists = set()
    for state in basis.stateList:
        dlists.update(gendlists(state, 3))

    for dlist in dlists:
        k1, k2, k3 = dlist

        k4 = k1+k2+k3
        clist = (k4,)

        V13.append((dlist,[clist]))


    V13 = Operator(V13, 3, 1, helper)


    V22 = []

    dlists = set()
    for state in basis.stateList:
        dlists.update(gendlists(state, 2))

    for dlist in dlists:
        (k1,k2) = dlist
        V22.append((dlist,[]))

        for k3 in range(max(-nmax+k1+k2,-nmax),
                min(int(floor((k1+k2)/2)),nmax)+1):

            k4 = k1+k2-k3
            clist = (k3,k4)

            # We always want to decrease the energy
            if helper.oscEnergy(clist) <= helper.oscEnergy(dlist) and\
                    helper.oscEnergy(clist) <= Emax+tol:

                V22[-1][1].append(clist)


    V22 = Operator(V22, 2, 2, helper)


    # NOTE V31 and V40 cannot increase the energy

    return V04, V13, V22
