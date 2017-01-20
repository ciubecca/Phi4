import gc
from sys import getsizeof as sizeof
import scipy
from math import factorial, floor, sqrt
import statefuncs
from statefuncs import Basis
from collections import Counter
import itertools
from statefuncs import Helper
from itertools import combinations, islice, permutations
from scipy import exp, pi, array
from scipy.special import binom
import bisect

tol = 10**(-10)


def filterDlist(dlist, nd, ntot, nmax):
    # TODO This can be sped up with the n-SUM algorithm
    if nd==ntot:
        return sum(dlist)==0
    elif nd==ntot-1:
        return abs(sum(dlist))<=nmax
    else:
        return True


# TODO The speed optimization resulting from using the N-Sum algorithm could be important
# @profile
def gendlists(state, nd, ntot, nmax):
    """ Generates a list of all the possible combinations of momenta in the state that
    can be annihilated
    state: input state in representation 1
    nd: number of annihilation operators (number of modes to annihilate)
    ntot: total number of annihilation and creation operators
    nmax: maximal wavenumber of the "lookup" basis
    """

    # XXX no call to list?
    x = itertools.chain.from_iterable(([n]*Zn for n,Zn in state))

    dlists = set(tuple(y) for y in combinations(x,nd))
    # XXX returns a generator expression
    return (dlist for dlist in dlists if filterDlist(dlist, nd, ntot, nmax))


def repr1torepr3(state):
    return list(itertools.chain.from_iterable([n]*Zn for n,Zn in state))



def pickN(state, N, occ):
    """ Pick all possible subset of particles of a given state """
    if N==0:
        return [[]]
    n,Zn = state[0]

    r = range(max(0,N-occ+Zn), min(Zn+1,N+1))
    return [[(n,Zm)]+x for Zm in r for x in pickN(state[1:],N-Zm,occ-Zn)]


def pickNM(state, N, M, occ):
    """ Divide a state into two subsets of particles """
    if state == []:
        return [[]]
    n,Zn = state[0]
    if N+M == Zn:
        return [[(n,N,M)]]

    r = range(max(0,N-occ+Zn), min(N+1,Zn+1))
    return [[(n,Zm,Zn-Zm)]+x for Zm in r for x in pickNM(state[1:],N-Zm,M-Zn+Zm,occ-Zn)]



def gendlistPairs(state, ndPair, ntotPair, nmax):

    ndTot = sum(ndPair)
    occn = statefuncs.occn(state)
    if occn < ndTot:
        return []

    ret = []

    for dlistTot in pickN(state, ndTot, occn):

        for dlist in pickNM(dlistTot, ndPair[0], ndPair[1], ndTot):

            dlist1 = tuple(itertools.chain.from_iterable([[n]*Zn1 for n,Zn1,Zn2 in dlist]))
            dlist2 = tuple(itertools.chain.from_iterable([[n]*Zn2 for n,Zn1,Zn2 in dlist]))

            if filterDlist(dlist1,ndPair[0],ntotPair[0],nmax) and\
                filterDlist(dlist2,ndPair[1],ntotPair[1],nmax):
                ret.append((dlist1, dlist2))

    return ret



def bose(x):
    """ computes the Bose factor of a product of oscillators  """
    return factorial(len(x))/scipy.prod(list(map(factorial,Counter(x).values())))


# XXX Check
parityFactors = [[1, sqrt(2)],[1/sqrt(2),1]]


class LocOperator():
    """
    Collection of oscillators with fixed number of creation and annihilation operators
    This is convenient to compute matrix elements and generate the high-energy basis
    from a set of tails
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

# This is just a dict where the keys are tuples of annihilation momenta and the
# values are indexes
        self.dlistPos = {}

# This is a list of lists. Each list contains all the possible creation momenta
# corresponding to a given set of annihilation momenta
        self.oscList = []

# List with same lenght as oscList. It contains the total energy of the corresponding
# oscillators
        self.oscEnergies = []

# Combinatorial and phase-space factors of the oscillators
        self.oscFactors = []


        for i, (dlist,clists) in enumerate(oscillators):
            clists = list(sorted(clists,key=helper.oscEnergy))

            self.dlistPos[dlist] = i

            self.oscList.append([self.torepr1(clist,dlist)
                for clist in clists])

            self.oscEnergies.append([helper.oscEnergy(clist)-helper.oscEnergy(dlist)
                for clist in clists])

            self.oscFactors.append([bose(clist)*bose(dlist)*binom(nc+nd,nc)\
                    *scipy.prod([1/sqrt(2*omega(n)*L) for n in clist+dlist])
                    for clist in clists])

    def torepr1(self, clist, dlist):
        """ This generates a list of tuples of the form [(n, Zc, Zd),...] from two separate
        tuples of the form (k1,...,kn) and (q1,...,qm), where the k's and q's are respectively
        the creation and annihilation momenta
        Zc and Zd are respectively the number of creation and annihilation operators at
        wavenumber n """

        wnlist = set(clist+dlist)
        cosc = Counter(clist)
        dosc = Counter(dlist)
        return [(n,cosc.get(n,0),dosc.get(n,0)) for n in wnlist]

    # @profile
    def computeMatrixElements(self, basis, i, lookupbasis, helper, statePos, Erange,
            ignKeyErr=False):
        """ Compute the matrix elements by applying all the oscillators in the operator
        to an element in the basis
        basis: set of states on which the operator acts
        i: index of the state in the basis
        lookupbasis: basis of states corresponding to the column indexes of the matrix
        helper: Helper instance
        statePos: dictionary where the keys are states in representation 2 (in tuple form)
        and the values are their position in the basis
        ignKeyErr: this must be set to True if the action of an oscillators on an input state
        can generate a state not in lookupbasis. This applies only in the computation of Vhh.
        Otherwise it should be set to False
        """

        cdef double x

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
        nmax = helper.nmax
        normFactors = helper.normFactors
        Emin, Emax = Erange

        # cycle over all the sets of momenta that can be annihilated
# XXX Check: we replaced lookupbasis.helper.nmax with helper.nmax
        for dlist in gendlists(state, self.nd, self.nd+self.nc, nmax):

            k = self.dlistPos[dlist]

# Only select the oscillators such that the sum of the state and oscillator energies
# lies within the bounds of the lookupbasis energies
            imin = bisect.bisect_left(self.oscEnergies[k], Emin-e-tol)
            imax = bisect.bisect_left(self.oscEnergies[k], Emax-e+tol)

            if imax <= imin:
                continue

            oscFactors = self.oscFactors[k][imin:imax]

            for i, osc in enumerate(self.oscList[k][imin:imax]):
                newstatevec = statevec[:]

                x = oscFactors[i]

                for n,Zc,Zd in osc:
                    newstatevec[n+nmax] += Zc-Zd
                    x *= normFactors[Zc, Zd, statevec[n+nmax]]

                if ignKeyErr:
                    try:
                        j = statePos[bytes(newstatevec)]
                    except KeyError:
                        continue
                else:
                    j = statePos[bytes(newstatevec)]

                x *= parityFactors[p][parityList[j]]
                data.append(x)
                col.append(j)

        return col, data


    # Generate high energy Hilbert space Hh from low energy Hilbert space Hl
    # as Hh = V*Hl
    def yieldBasis(self, basis, EL):
        """ Return a set of tuples of representation 2 states, all of which are not
        connected by spatial parity transformations.
        basis: set of "selected" low energy states on which to act
        EL: maximal energy of the generated high-energy states
        """

        nmax = self.helper.nmax

        for i, state in enumerate(basis):

            statevec = self.helper.torepr2(state)
            e = basis.energyList[i]

            for dlist in gendlists(state, self.nd, self.nd+self.nc, nmax):
                k = self.dlistPos[dlist]

                imax = bisect.bisect_left(self.oscEnergies[k], EL-e+tol)

                for i, osc in enumerate(self.oscList[k][:imax]):
                    newstatevec = statevec[:]
                    for n,Zc,Zd in osc:
                        newstatevec[n+nmax] += Zc-Zd
                    yield bytes(newstatevec)

        # print(sizeof(tuple(statevec)))
        # print(sizeof(tuple(statevec)[1]))
        # print("type of data in statevec:", type(statevec[0]))
        # print("type of data in tuples", type(t1[0]))



class BilocOperator():
    def __init__(self, JointOscList, ndPair, ncPair, helper):
        """
        oscillators: list of tuples. The first element of the tuple is a tuple of
        wavenumbers of annihilation operators, and the second element a list of
        tuples of wavenumbers of creation operators
        basis: basis on which the operator will act
        nd: number of annihilation operators
        nc: number of creation operators
        """

        self.ndPair = ndPair
        self.ncPair = ncPair
        self.ntotPair = tuple((nd+nc for nc,nd in zip(ncPair,ndPair)))
        omega = helper.omega
        L = helper.L
        m = helper.m
        self.helper = helper

        self.dlistPairPos = {}
        self.oscList = []
        self.oscEnergies = []
        self.oscFactors = []

        oscEnergy = helper.oscEnergy

        def pairOscEnergy(pair):
            return oscEnergy(pair[0])+oscEnergy(pair[1])

        for i, (dlistPair,clistPairs) in enumerate(JointOscList):
            clistPairs = list(sorted(clistPairs, key=pairOscEnergy))

            self.dlistPairPos[dlistPair] = i

            self.oscList.append([self.torepr1(clistPair,dlistPair)
                for clistPair in clistPairs])

            self.oscEnergies.append([pairOscEnergy(clistPair)-pairOscEnergy(dlistPair)
                for clistPair in clistPairs])

            self.oscFactors.append([
                scipy.prod([bose(clist)*bose(dlist)*binom(nc+nd,nc)*\
                    scipy.prod([1/sqrt(2*omega(n)*L) for n in clist+dlist])
                    for clist,dlist,nc,nd in zip(clistPair,dlistPair,ncPair,ndPair)])
                        for clistPair in clistPairs])


    def torepr1(self, clistPair, dlistPair):
        """ This generates a list of tuples of the form [(n, Zc, Zd),...] from two separate
        tuples of the form (k1,...,kn) and (q1,...,qm), where the k's and q's are respectively
        the creation and annihilation momenta
        Zc and Zd are respectively the number of creation and annihilation operators at
        wavenumber n """
        clist = tuple(sorted(clistPair[0]+clistPair[1]))
        dlist = tuple(sorted(dlistPair[0]+dlistPair[1]))
        wnlist = set(clist+dlist)
        cosc = Counter(clist)
        dosc = Counter(dlist)
        return list((n,cosc.get(n,0),dosc.get(n,0)) for n in wnlist)


    def computeMatrixElements(self, basis, i, lookupbasis, helper, statePos, Erange,
                                ignKeyErr=True):
        """ Compute the matrix elements by applying all the oscillators in the operator
        to an element in the basis
        basis: set of states on which the operator acts
        i: index of the state in the basis
        lookupbasis: basis of states corresponding to the column indexes of the matrix
        helper: Helper instance
        statePos: dictionary where the keys are states in representation 2 (in tuple form)
        and the values are their position in the basis
        ignKeyErr: this must be set to True if the action of an oscillators on an input state
        can generate a state not in lookupbasis. This applies only in the computation of Vhh.
        Otherwise it should be set to False
        """


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
        nmax = helper.nmax
        Emin, Emax = Erange

        normFactors = helper.normFactors

        # cycle over all the sets of momenta that can be annihilated
        for dlistPair in gendlistPairs(state, self.ndPair, self.ntotPair,
                lookupbasis.helper.nmax):

            k = self.dlistPairPos[dlistPair]

# Only select the oscillators such that the sum of the state and oscillator energies
# lies within the bounds of the lookupbasis energies
            imin = bisect.bisect_left(self.oscEnergies[k], Emin-e-tol)
            imax = bisect.bisect_left(self.oscEnergies[k], Emax-e+tol)

            if imax <= imin:
                continue

            oscFactors = self.oscFactors[k][imin:imax]

            for i, osc in enumerate(self.oscList[k][imin:imax]):
                newstatevec = statevec[:]

                x = oscFactors[i]

                for n,Zc,Zd in osc:
                    newstatevec[n+nmax] += Zc-Zd
                    x *= normFactors[Zc, Zd, statevec[n+nmax]]

                if ignKeyErr:
                    try:
                        j = statePos[tuple(newstatevec)]
                    except KeyError:
                        continue
                else:
                    j = statePos[tuple(newstatevec)]

                x *= parityFactors[p][parityList[j]]
                data.append(x)
                col.append(j)


        return col, data



def V4OpsHalf(basis):
    """ Generate half of the oscillators () of the V4 operator
    basis: basis of all the low-energy states below ET """

    nmax = basis.nmax
    Emax = basis.Emax
    helper = basis.helper

    dlist = ()
# The list of annihilation momenta is empty
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

# Generate an LocOperator instance from the computed set of oscillators
    V40 = LocOperator(V40, 0, 4, helper)


    V31 = []
    for k1 in range(-nmax,nmax+1):
# The set of annihilation momenta contains just one momentum
        dlist = (k1,)
        V31.append((dlist,[]))

        for k2 in range(-nmax,nmax+1):
            for k3 in range(max(-nmax+k1-k2,k2),
                           min(int(floor((k1-k2)/2)),nmax)+1):

                k4 = k1-k2-k3
                clist = (k2,k3,k4)

                if helper.oscEnergy(clist) <= Emax+tol:
                    V31[-1][1].append(clist)

    V31 = LocOperator(V31, 1, 3, helper)


    V22 = []
    for k1 in range(-nmax,nmax+1):
        for k2 in range(k1,nmax+1):
# The set of annihilation momenta is a pair
            dlist = (k1,k2)
            V22.append((dlist,[]))

            for k3 in range(max(-nmax+k1+k2,-nmax),
                    min(int(floor((k1+k2)/2)),nmax)+1):

                k4 = k1+k2-k3
                clist = (k3,k4)

                if helper.oscEnergy(clist) <= Emax+tol and\
                    sorted([abs(k3),abs(k4)])<=sorted([abs(k1),abs(k2)]):
                    # only consider lexicographically ordered part of V22
                    # but also including diagonal part which will be separated below

                    V22[-1][1].append(clist)


    V22 = LocOperator(V22, 2, 2, helper)

    return V40, V31, V22



def V2OpsHalf(basis):
    """ Generate half of the oscillators () of the V2 operator
    basis: basis of all the low-energy states below ET """

    nmax = basis.nmax
    Emax = basis.Emax
    helper = basis.helper

    dlist = ()
    V20 = [(dlist, [])]
    for k1 in range(-nmax,1):
        k2 = -k1
        clist = (k1,k2)

        V20[-1][1].append(clist)

    V20 = LocOperator(V20, 0, 2, helper)


    V11 = []
    for k1 in range(-nmax,nmax+1):
        dlist = (k1,)
        V11.append((dlist,[]))

        k2 = k1
        clist = (k2,)
        V11[-1][1].append(clist)

    V11 = LocOperator(V11, 1, 1, helper)

    return V20, V11




def V4OpsSelectedHalf(basis, Emax, idxList=None):
    """ Selected set of oscillators of half of the V4 operator between selected states
    basis: basis which is acted upon
    Emin: minimum energy of the states to be generated
    Emax: maximal energy of states to be generated
    idxList: subset of indices of the basis which is acted upon
    """

    helper = Helper(basis.helper.m, basis.helper.L, max(Emax,basis.Emax))
    nmax = helper.nmax
    oscEnergy = helper.oscEnergy

    if idxList == None:
        idxList = range(basis.size)

    opsList = []

    for nd in (0,1,2):
        nc = 4-nd

        dlists = gendlistsfromBasis(basis, idxList, nmax, nd, 4)
        oscList = []

        for dlist in dlists:
            clists = [clist for clist in createClistsV4(nmax, dlist, nc) if
                    oscEnergy(clist) <= Emax+tol]
            if nd==2:
                clists = [clist for clist in clists if
                    sorted([abs(dlist[0]),abs(dlist[1])])<=
                    sorted([abs(clist[0]),abs(clist[1])]) ]

            oscList.append((dlist, clists))

        opsList.append(LocOperator(oscList,nd,nc,helper=helper))

    return opsList


def V2OpsSelectedFull(basis, Emax):
    """ Selected set of oscillators between some selected low-energy states
    and states with energy <= Emax
    """

    helper = Helper(basis.helper.m, basis.helper.L, max(Emax,basis.Emax))
    nmax = helper.nmax
    oscEnergy = helper.oscEnergy

    opsList = []

    for nd in (0,1,2):
        nc = 2-nd

        dlists = gendlistsfromBasis(basis, range(basis.size), nmax, nd, 2)
        oscList = []

        for dlist in dlists:
            clists = [clist for clist in createClistsV2(nmax, dlist, nc) if
                    oscEnergy(clist) <= Emax+tol]
            oscList.append((dlist, clists))

        opsList.append(LocOperator(oscList,nd,nc,helper=helper))

    return opsList


def V4OpsSelectedFull(basis, Emax, idxList=None):
    """ Selected set of oscillators of the full V4 operator between some selected states
    basis: basis which is acted upon
    Emin: minimum energy of the states to be generated
    Emax: maximal energy of states to be generated
    idxList: subset of indices of the basis which is acted upon
    """

    helper = Helper(basis.helper.m, basis.helper.L, max(Emax,basis.Emax))
    nmax = helper.nmax
    oscEnergy = helper.oscEnergy

    if idxList == None:
        idxList = range(basis.size)

    opsList = []

    for nd in (0,1,2,3,4):
        nc = 4-nd

        dlists = gendlistsfromBasis(basis, idxList, nmax, nd, 4)
        oscList = []

        for dlist in dlists:
            clists = [clist for clist in createClistsV4(nmax, dlist, nc) if
                    oscEnergy(clist) <= Emax+tol]
            oscList.append((dlist, clists))

        opsList.append(LocOperator(oscList,nd,nc,helper=helper))

    return opsList


def V6OpsSelectedFull(basis, Emax):
    """ Selected set of oscillators between some selected low-energy states
    and states with energy <= Emax
    """

    helper = Helper(basis.helper.m, basis.helper.L, max(Emax,basis.Emax))
    nmax = helper.nmax
    oscEnergy = helper.oscEnergy

    opsList = []

    for nd in (0,1,2,3,4,5,6):
        nc = 6-nd

        dlists = gendlistsfromBasis(basis, range(basis.size), nmax, nd, 6)
        # print(dlists)
        oscList = []

        for dlist in dlists:
            clists = [clist for clist in createClistsV6(nmax, dlist, nc) if
                    oscEnergy(clist) <= Emax+tol]
            oscList.append((dlist, clists))

        opsList.append(LocOperator(oscList,nd,nc,helper=helper))

    return opsList



def gendlistsfromBasis(basis, idxList, nmax, nd, ntot):
    ret = set()

    for i in idxList:
        state = basis.stateList[i]
        ret.update(gendlists(state=state, nd=nd, ntot=ntot, nmax=nmax))
    return ret



def createClistsV2(nmax, dlist, nc):

    if len(dlist) != 2-nc:
        raise ValueError
    if nc==0:
        return [()]
    elif nc==1:
        return [(sum(dlist),)]
    elif nc==2:
        return [(k3,-k3) for k3 in range(-nmax,1)]



def createClistsV4(nmax, dlist, nc):

    if len(dlist) != 4-nc:
        raise ValueError
    clists = []

    if nc==0:
        clists.append(())
    elif nc==1:
        clists.append((sum(dlist),))
    elif nc==2:
        k1,k2 = dlist
        for k3 in range(max(-nmax+k1+k2,-nmax),
                min(int(floor((k1+k2)/2)),nmax)+1):

            k4 = k1+k2-k3
            clists.append((k3,k4))

    elif nc==3:
        (k1,) = dlist
        for k2 in range(-nmax,nmax+1):
            for k3 in range(max(-nmax+k1-k2,k2),
                           min(int(floor((k1-k2)/2)),nmax)+1):

                k4 = k1-k2-k3
                clists.append((k2,k3,k4))

    elif nc==4:
        clists = []
        for k1 in range(-nmax,nmax+1):
            for k2 in range(k1,nmax+1):
                # NOTE the boundaries for k3 ensure that k3<=k4<=nmax
                for k3 in range(max(-nmax-k1-k2,k2),
                        min(int(floor((-k1-k2)/2)),nmax)+1):
                    k4 = -k1-k2-k3
                    clists.append((k1,k2,k3,k4))

    return clists



def createClistsV6(nmax, dlist, nc):

    if len(dlist) != 6-nc:
        raise ValueError
    clists = []

    if nc==0:
        clists.append(())
    elif nc==1:
        clists.append((sum(dlist),))

    elif nc==2:
        sumk = sum(dlist)
        for k5 in range(-nmax, nmax+1):
            k6 = sumk-k5
            if k5<=k6<=nmax:
                clists.append((k5,k6))

    elif nc==3:
        sumk = sum(dlist)
        for k4 in range(-nmax,nmax+1):
            for k5 in range(k4,nmax+1):
                k6 = sumk-k4-k5
                if k5<=k6<=nmax:
                    clists.append((k4,k5,k6))

    elif nc==4:
        sumk = sum(dlist)
        for k3 in range(-nmax,nmax+1):
            for k4 in range(k3,nmax+1):
                for k5 in range(k4,nmax+1):
                    k6 = sumk-k3-k4-k5
                    if k5<=k6<=nmax:
                        clists.append((k3,k4,k5,k6))
    elif nc==5:
        sumk = sum(dlist)
        for k2 in range(-nmax,nmax+1):
            for k3 in range(k2,nmax+1):
                for k4 in range(k3,nmax+1):
                    for k5 in range(k4,nmax+1):
                        k6 = sumk-k2-k3-k4-k5
                        if k5<=k6<=nmax:
                            clists.append((k2,k3,k4,k5,k6))
    elif nc==6:
        for k1 in range(-nmax,nmax+1):
            for k2 in range(k1,nmax+1):
                for k3 in range(k2,nmax+1):
                    for k4 in range(k3,nmax+1):
                        for k5 in range(k4,nmax+1):
                            k6 = -k1-k2-k3-k4-k5
                            if k5<=k6<=nmax:
                                clists.append((k1,k2,k3,k4,k5,k6))
    return clists



def gendlistPairsfromBasis(basis, nmax, ndPair, ntotPair):
    ret = set()

    # print("nmax", nmax)
    # print("ndPair", ndPair)
    # print("ntotPair", ntotPair)

    for state in basis:
        # print(state)
        ret.update(gendlistPairs(state=state, ndPair=ndPair,
            ntotPair=ntotPair, nmax=nmax))
    return ret



def V2V4Ops(basis):
    """
    Bilocal operators of :V2 V4:
    """

    helper = basis.helper
    nmax = helper.nmax
    Emax = basis.Emax
    oscEnergy = helper.oscEnergy

    ntotPair = (2,4)

    createClistsV = {2:createClistsV2, 4:createClistsV4}

    opsList = []

    for nd1 in (0,1,2):
        for nd2 in (0,1,2,3,4):
            ndPair = (nd1,nd2)

            ncPair = tuple(ntot-nd for ntot,nd in zip(ntotPair,ndPair))

            # print(nd1, nd2)
            dlistPairs = gendlistPairsfromBasis(basis, nmax, ndPair, ntotPair)

            JointOscList = []

            for dlistPair in dlistPairs:

                x1 = createClistsV[2](nmax, dlistPair[0], ncPair[0])
                x2 = createClistsV[4](nmax, dlistPair[1], ncPair[1])


                clistPairs = [(clist1,clist2) for clist1 in x1 for clist2 in x2
                        if oscEnergy(clist1)+oscEnergy(clist2) <= Emax+tol]

                JointOscList.append((dlistPair, clistPairs))

            opsList.append(BilocOperator(JointOscList,ndPair,ncPair,helper=helper))

    return opsList




def V4V4Ops(basis):
    """
    Bilocal operators of :V4 V4:
    """

    helper = basis.helper
    nmax = helper.nmax
    Emax = basis.Emax
    oscEnergy = helper.oscEnergy

    ntotPair = (4,4)

    opsList = []

    for nd1 in (0,1,2,3,4):
        for nd2 in (0,1,2,3,4):
            ndPair = (nd1,nd2)

            ncPair = tuple(ntot-nd for ntot,nd in zip(ntotPair,ndPair))

            dlistPairs = gendlistPairsfromBasis(basis, nmax, ndPair, ntotPair)

            JointOscList = []

            for dlistPair in dlistPairs:

                x1 = createClistsV4(nmax, dlistPair[0], ncPair[0])
                x2 = createClistsV4(nmax, dlistPair[1], ncPair[1])


                clistPairs = [(clist1,clist2) for clist1 in x1 for clist2 in x2
                        if oscEnergy(clist1)+oscEnergy(clist2) <= Emax+tol]

                JointOscList.append((dlistPair, clistPairs))

            opsList.append(BilocOperator(JointOscList,ndPair,ncPair,helper=helper))

    return opsList

