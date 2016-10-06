import scipy
from math import factorial, floor, sqrt
from statefuncs import Basis
from collections import Counter
import itertools
from statefuncs import Helper
from itertools import combinations, islice
from scipy import exp, pi, array
import bisect

tol = 10**(-10)

def gendlists(state, nd, ntot, nmax):
    x = list(itertools.chain.from_iterable([[n]*Zn for n,Zn in state]))
    ret = set(tuple(y) for y in combinations(x,nd))
    # TODO This can be sped up with the n-SUM algorithm
    if nd==ntot:
        return set(dlist for dlist in ret if sum(dlist)==0)
    else if nd==ntot-1:
        return set(dlist for dlist in ret if abs(sum(dlist))<=nmax)
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


    def computeMatrixElements(self, basis, i, lookupbasis, statePos, helper,
            ignoreKeyError=False):
        # List of columns indices of generated basis elements
        col = []
        # List of partial matrix elements
        data = []

        # I define these local variables outside the loops for performance reasons
        e = basis.energyList[i]
        p = basis.parityList[i]
        state = basis.stateList[i]

        statevec = helper.torepr2(state)
        parityList = lookupbasis.parityList
        Emin = lookupbasis.Emin
        Emax = lookupbasis.Emax
        nmax = helper.nmax

        normFactors = helper.normFactors

        # cycle over all the sets of momenta that can be annihilated
        for dlist in gendlists(state, self.nd):

            try:
                k = self.dlistPos[dlist]
            except KeyError as err:
                print(helper.oscEnergy(dlist))
                raise err

            imin = bisect.bisect_left(self.oscEnergies[k], Emin-e-tol)
            imax = bisect.bisect_left(self.oscEnergies[k], Emax-e+tol)

            if imax <= imin:
                continue

            # datapart = (self.oscFactors[k][imin:imax])[:]
            oscFactors = self.oscFactors[k][imin:imax]

            # This cycle could be moved to C
            for i, osc in enumerate(self.oscList[k][imin:imax]):
                newstatevec = statevec[:]

                x = oscFactors[i]

                for n,Zc,Zd in osc:
                    newstatevec[n+nmax] += Zc-Zd
                    # datapart[i] *= normFactors[Zc, Zd, statevec[n+nmax]]
                    x *= normFactors[Zc, Zd, statevec[n+nmax]]

                # This is in case the lookup basis is "selected"
                if ignoreKeyError:
                    try:
                        j = statePos[tuple(newstatevec)]
                    except KeyError:
                        continue
                else:
                    j = statePos[tuple(newstatevec)]

                x *= parityFactors[p][parityList[j]]
                data.append(x)
                col.append(j)

                # datapart[i] *= parityFactors[p][parityList[j]]
                # col.append(j)

            # data += datapart

        return col, data


    # Generate high energy Hilbert space Hh from low energy Hilbert space Hl
    # as Hh = V*Hl
    def genBasis(self, basis, ET, EL):

        nmax = self.helper.nmax
        stateset = set()

        for i, state in enumerate(basis):

            statevec = self.helper.torepr2(state)
            e = basis.energyList[i]

            for dlist in gendlists(state, self.nd):
                k = self.dlistPos[dlist]

                imin = bisect.bisect_left(self.oscEnergies[k], ET-e+tol)
                imax = bisect.bisect_left(self.oscEnergies[k], EL-e+tol)
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



def V4Operators(basis):

    nmax = basis.nmax
    Emax = basis.Emax
    helper = basis.helper

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



def V2Operators(basis):
    """
    basis: set of states on which the operators act
    """
    nmax = basis.nmax
    Emax = basis.Emax
    helper = basis.helper

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


def V4Operatorshl(basis, EL):
    """
    basis: set of states on which the operators act
    EL: maximal energy of the states corresponding to column indices of the matrix
    """
    helper = Helper(basis.helper.m, basis.helper.L, EL)
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

                if helper.oscEnergy(clist) <= EL+tol:
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

                if helper.oscEnergy(clist) <= EL+tol\
                    and helper.oscEnergy(clist) > helper.oscEnergy(dlist):
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
            if helper.oscEnergy(clist) > helper.oscEnergy(dlist) and\
                    helper.oscEnergy(clist) <= EL+tol:

                V22[-1][1].append(clist)


    V22 = Operator(V22, 2, 2, helper)


    # NOTE V13 and V04 cannot increase the energy

    return V40, V31, V22


# Generate all the oscillators between the "selected" high-energy basis and the
# low-energy basis
def V4OperatorsLh(basis, Emax):
    """
    basis: set of states on which the operators act
    Emax: maximal energy of the states corresponding to columns indices of the matrix
    """

    # These two lines of code are just used to compute nmax of the creation operators
    helper = Helper(basis.helper.m, basis.helper.L, Emax)
    nmax = helper.nmax

    helper = basis.helper

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

        # This can be improved
        if -nmax <= k4 <= nmax:
            V13.append((dlist,[clist]))
        else:
            V13.append((dlist,[]))

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


def V4Operatorshh(basis):

    helper = basis.helper
    nmax = helper.nmax
    Emin = basis.Emin
    Emax = basis.Emax

    dlist = ()
    V40 = [(dlist, [])]
    for k1 in range(-nmax,nmax+1):
        for k2 in range(k1,nmax+1):
            # NOTE the boundaries for k3 ensure that k3<=k4<=nmax
            for k3 in range(max(-nmax-k1-k2,k2),
                    min(int(floor((-k1-k2)/2)),nmax)+1):

                k4 = -k1-k2-k3
                clist = (k1,k2,k3,k4)

                if helper.oscEnergy(clist) <= Emax-Emin+tol:
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

                if helper.oscEnergy(clist) <= Emax+tol:
                    V31[-1][1].append(clist)

    V31 = Operator(V31, 1, 3, helper)


    V22 = []

    dlists = set()
    for state in basis.stateList:
        dlists.update(gendlists(state, 2))

    for dlist in dlists:
        k1,k2 = dlist
        V22.append((dlist,[]))

        for k3 in range(max(-nmax+k1+k2,-nmax),
                min(int(floor((k1+k2)/2)),nmax)+1):

            k4 = k1+k2-k3
            clist = (k3,k4)

            # NOTE The check on dlist is useless here but it's needed
            # if we generalize the code to other operators
            if helper.oscEnergy(clist) <= Emax+tol and\
                sorted([abs(k3),abs(k4)])<=sorted([abs(k1),abs(k2)]):
                # only consider lexicographically ordered part of V22
                # but also including diagonal part which will be separated below

                V22[-1][1].append(clist)


    V22 = Operator(V22, 2, 2, helper)

    return V40, V31, V22




def V6Operatorshh(basis):

    helper = basis.helper
    nmax = helper.nmax
    Emin = basis.Emin
    Emax = basis.Emax

    dlist = ()
    V40 = [(dlist, [])]
    for k1 in range(-nmax,nmax+1):
        for k2 in range(k1,nmax+1):
            # NOTE the boundaries for k3 ensure that k3<=k4<=nmax
            for k3 in range(max(-nmax-k1-k2,k2),
                    min(int(floor((-k1-k2)/2)),nmax)+1):

                k4 = -k1-k2-k3
                clist = (k1,k2,k3,k4)

                if helper.oscEnergy(clist) <= Emax-Emin+tol:
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

                if helper.oscEnergy(clist) <= Emax+tol:
                    V31[-1][1].append(clist)

    V31 = Operator(V31, 1, 3, helper)


    V22 = []

    dlists = set()
    for state in basis.stateList:
        dlists.update(gendlists(state, 2))

    for dlist in dlists:
        k1,k2 = dlist
        V22.append((dlist,[]))

        for k3 in range(max(-nmax+k1+k2,-nmax),
                min(int(floor((k1+k2)/2)),nmax)+1):

            k4 = k1+k2-k3
            clist = (k3,k4)

            # NOTE The check on dlist is useless here but it's needed
            # if we generalize the code to other operators
            if helper.oscEnergy(clist) <= Emax+tol and\
                sorted([abs(k3),abs(k4)])<=sorted([abs(k1),abs(k2)]):
                # only consider lexicographically ordered part of V22
                # but also including diagonal part which will be separated below

                V22[-1][1].append(clist)


    V22 = Operator(V22, 2, 2, helper)

    return V40, V31, V22

