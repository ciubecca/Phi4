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

# @profile
# TODO The speed optimization resulting from using the N-Sum algorithm could be important
def gendlists(state, nd, ntot, nmax):
    """ Generates a list of all the possible combinations of momenta in the state that
    can be annihilated
    state: input state in representation 1
    nd: number of annihilation operators (number of modes to annihilate)
    ntot: total number of annihilation and creation operators
    nmax: maximal wavenumber of the "lookup" basis
    """

    x = list(itertools.chain.from_iterable([[n]*Zn for n,Zn in state]))

    # TODO This can be sped up with the n-SUM algorithm
    if nd==ntot:
        # In this case the momenta must sum to 0
        ret = set(tuple(y) for y in combinations(x,nd))
        return set(dlist for dlist in ret if sum(dlist)==0)
    elif nd==ntot-1:
# In this case we can already exclude momenta whose sum exceeds nmax
        ret = set(tuple(y) for y in combinations(x,nd))
        return set(dlist for dlist in ret if abs(sum(dlist))<=nmax)
    else:
        ret = set(tuple(y) for y in combinations(x,nd))
        return ret

def bose(x):
    """ computes the Bose factor of a product of oscillators  """
    return factorial(len(x))/scipy.prod(list(map(factorial,Counter(x).values())))


# XXX Check
parityFactors = [[1, sqrt(2)],[1/sqrt(2),1]]

class LocOperator(Operator):

class BilocOperator(Operator):

class Operator():
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

            self.oscList.append([self.torepr1(clist,dlist) for clist in clists])

            self.oscEnergies.append([helper.oscEnergy(clist)-helper.oscEnergy(dlist)
                for clist in clists])

            self.oscFactors.append([bose(clist)*bose(dlist)*scipy.special.binom(nc+nd,nc)\
                    *scipy.prod([1/sqrt(2*omega(n)*L) for n in clist+dlist])
                    for clist in clists])


    def torepr1(clist, dlist):
        """ This generates a list of tuples of the form [(n, Zc, Zd),...] from two separate
        tuples of the form (k1,...,kn) and (q1,...,qm), where the k's and q's are respectively
        the creation and annihilation momenta
        Zc and Zd are respectively the number of creation and annihilation operators at
        wavenumber n """

        wnlist = set(clist+dlist)
        cosc = Counter(clist)
        dosc = Counter(dlist)
        return list((n,cosc.get(n,0),dosc.get(n,0)) for n in wnlist)


    # @profile
    def computeMatrixElements(self, basis, i, lookupbasis, helper, statePos,
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
        for dlist in gendlists(state, self.nd, self.nd+self.nc, lookupbasis.helper.nmax):

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
                        j = statePos[tuple(newstatevec)]
                    except KeyError:
                        continue
                else:
                    j = statePos[tuple(newstatevec)]

                x *= parityFactors[p][parityList[j]]
                data.append(x)
                col.append(j)


        return col, data


    # Generate high energy Hilbert space Hh from low energy Hilbert space Hl
    # as Hh = V*Hl
    def genBasis(self, basis, EL):
        """ Return a set of tuples of representation 2 states, all of which are not
        connected by spatial parity transformations.
        basis: set of "selected" low energy states on which to act
        EL: maximal energy of the generated high-energy states
        """

        nmax = self.helper.nmax
        stateset = set()

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
                    t1 = tuple(newstatevec)
                    t2 = tuple(newstatevec[::-1])
                    if (t1 not in stateset) and (t2 not in stateset):
                        stateset.add(t1)

        return stateset



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

# Generate an Operator instance from the computed set of oscillators
    V40 = Operator(V40, 0, 4, helper)


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

    V31 = Operator(V31, 1, 3, helper)


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


    V22 = Operator(V22, 2, 2, helper)

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

    V20 = Operator(V20, 0, 2, helper)


    V11 = []
    for k1 in range(-nmax,nmax+1):
        dlist = (k1,)
        V11.append((dlist,[]))

        k2 = k1
        clist = (k2,)
        V11[-1][1].append(clist)

    V11 = Operator(V11, 1, 1, helper)

    return V20, V11


# Oscillators between the selected states in basis and all the states in the range [0,Emax]
def V4OpsSelectedFull(basis, Emax):
    """ Returns Operator instances containing all the oscillators of the V4 operator
    acting on basis.  All the parts of the operator are computed.
    The lists of annihilation momenta are not all the possible ones, but just those
    extracted from the particular states in the basis.
    This is method is used to compute Vlh and VhL
    basis: "selected" basis of states (e.g. low-energy tails)
    Emax: maximal energy of the states generated via the operator.
    """

    helper = Helper(basis.helper.m, basis.helper.L, max(Emax,basis.Emax))
    nmax = helper.nmax

    #############
    # a^ a^ a^ a^
    #############
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

    #############
    # a^ a^ a^ a
    #############
    V31 = []

    dlists = set()
# Generate all the lists of annihilation momenta from the states in the basis
    for state in basis.stateList:
        dlists.update(gendlists(state, 1, 4, nmax))

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

    #############
    # a^ a^ a a
    #############
    V22 = []

    dlists = set()
    for state in basis.stateList:
        dlists.update(gendlists(state, 2, 4, nmax))

    for dlist in dlists:
        (k1,k2) = dlist
        V22.append((dlist,[]))

        for k3 in range(max(-nmax+k1+k2,-nmax),
                min(int(floor((k1+k2)/2)),nmax)+1):

            k4 = k1+k2-k3
            clist = (k3,k4)

            if helper.oscEnergy(clist) <= Emax+tol:
                V22[-1][1].append(clist)

    V22 = Operator(V22, 2, 2, helper)


    #############
    # a^ a a a
    #############
    V13 = []

    dlists = set()
    for state in basis.stateList:
        dlists.update(gendlists(state, 3, 4, nmax))

    for dlist in dlists:
        k1, k2, k3 = dlist

        k4 = k1+k2+k3
        clist = (k4,)

        V13.append((dlist,[clist]))

    V13 = Operator(V13, 3, 1, helper)


    ############
    # a a a a
    ###########
    V04 = []

    dlists = set()
    for state in basis.stateList:
    # This has to be implemented differently
        dlists.update(gendlists(state, 4, 4, nmax))

    for dlist in dlists:
        V04.append((dlist,[()]))

    V04 = Operator(V04, 4, 0, helper)


    return V40, V31, V22, V13, V04



def V4OpsSelectedHalf(basis):
    """ Half of the operators of V4 acting on the "selected" basis of states
    (e.g. the set of "selected" high-energy tails). It takes into account only the
    annihilation part.
    Not all the possible annihilation momenta are included. This method is used
    to compute Vhh
    """

    helper = basis.helper
    nmax = helper.nmax
    Emax = basis.Emax
    oscEnergy = helper.oscEnergy

    V04 = []

    dlists = set()
    for state in basis.stateList:
        dlists.update(gendlists(state, 4, 4, nmax))

    for dlist in dlists:
        # XXX Some of these elements can be excluded
        V04.append((dlist,[()]))

    V04 = Operator(V04, 4, 0, helper)


    V13 = []

    dlists = set()
    for state in basis.stateList:
        dlists.update(gendlists(state, 3, 4, nmax))

    for dlist in dlists:
        k1, k2, k3 = dlist

        k4 = k1+k2+k3
        clist = (k4,)

        # XXX Some of these elements can be excluded
        V13.append((dlist,[clist]))

    V13 = Operator(V13, 3, 1, helper)


    V22 = []

    dlists = set()
    for state in basis.stateList:
        dlists.update(gendlists(state, 2, 4, nmax))

    for dlist in dlists:
        (k1,k2) = dlist
        V22.append((dlist,[]))

        for k3 in range(max(-nmax+k1+k2,-nmax),
                min(int(floor((k1+k2)/2)),nmax)+1):

            k4 = k1+k2-k3
            clist = (k3,k4)

            if oscEnergy(clist) <= Emax + tol \
                and sorted([abs(k3),abs(k4)])<=sorted([abs(k1),abs(k2)]):
                V22[-1][1].append(clist)


    V22 = Operator(V22, 2, 2, helper)

    return V04, V13, V22




def V6OpsSelectedHalf(basis):
    """ Half of the operators of V6 acting on the "selected" basis of states
    (e.g. the set of "selected" low-energy tails). It takes into account only the
    annihilation part.
    Not all the possible annihilation momenta are included. This method is used
    when computing the local V6 matrix on the low-energy tails
    """

    helper = basis.helper
    nmax = helper.nmax
    Emax = basis.Emax
    oscEnergy = helper.oscEnergy


    V06 = []

    dlists = set()
    for state in basis.stateList:
        dlists.update(gendlists(state, 6, 6, nmax))

    for dlist in dlists:
        # TODO Some of these elements can be excluded
        V06.append((dlist,[()]))

    V06 = Operator(V06, 6, 0, helper)


    V15 = []

    dlists = set()
    for state in basis.stateList:
        dlists.update(gendlists(state, 5, 6, nmax))

    for dlist in dlists:
        k1, k2, k3, k4, k5 = dlist

        k6 = k1+k2+k3+k4+k5
        clist = (k6,)

        # TODO Some of these elements can be excluded
        V13.append((dlist,[clist]))

    V15 = Operator(V15, 5, 1, helper)


    V24 = []

    dlists = set()
    for state in basis.stateList:
        dlists.update(gendlists(state, 2, 6, nmax))

    for dlist in dlists:
        (k1,k2,k3,k4) = dlist
        V24.append((dlist,[]))

        for k5 in range(max(-nmax+k1+k2+k3+k4,-nmax),
                min(int(floor((k1+k2+k3+k4)/2)),nmax)+1):

            k6 = k1+k2+k3+k4-k5
            clist = (k5,k6)

            if oscEnergy(clist) <= Emax+ tol:
                V22[-1][1].append(clist)


    V24 = Operator(V24, 4, 2, helper)



    V33 = []

    dlists = set()
    for state in basis.stateList:
        dlists.update(gendlists(state, 3, 6, nmax))

    for dlist in dlists:
        (k1,k2,k3) = dlist
        V33.append((dlist,[]))

        for k4 in range(-nmax, nmax+1):

            for k5 in range(max(-nmax+k1+k2+k3-k4,-nmax),
                min(int(floor((k1+k2+k3-k4)/2)),nmax)+1):

                k6 = k1+k2+k3-k4-k5
                clist = (k4,k5,k6)

                if oscEnergy(clist) <= Emax+ tol \
                    and sorted([abs(min(dlist)),abs(max(dlist))]) <= \
                        sorted([abs(min(clist)),abs(max(clist))]):
                            # XXX This is a generalization of the lexicographical
#sorting condition
                    V33[-1][1].append(clist)


    V33 = Operator(V33, 3, 3, helper)

    return V06, V15, V24, V33



def V2V4OpsSelectedHalf(basis):
    """ Half of the bilocal operators of :V2 V4: acting on the "selected" basis of states
    (e.g. the set of "selected" low-energy tails). It takes into account only the
    annihilation part.
    Not all the possible annihilation momenta are included. This method is used
    when computing the local V6 matrix on the low-energy tails
    """

    helper = basis.helper
    nmax = helper.nmax
    Emax = basis.Emax
    oscEnergy = helper.oscEnergy


    V06 = []

    dlists = set()
    for state in basis.stateList:
        dlists.update(gendlists(state, 6, 6, nmax))

    for dlist in dlists:
        # TODO Some of these elements can be excluded
        V06.append((dlist,[()]))

    V06 = Operator(V06, 6, 0, helper)


    V15 = []

    dlists = set()
    for state in basis.stateList:
        dlists.update(gendlists(state, 5, 6, nmax))

    for dlist in dlists:
        k1, k2, k3, k4, k5 = dlist

        k6 = k1+k2+k3+k4+k5
        clist = (k6,)

        # TODO Some of these elements can be excluded
        V13.append((dlist,[clist]))

    V15 = Operator(V15, 5, 1, helper)


    V24 = []

    dlists = set()
    for state in basis.stateList:
        dlists.update(gendlists(state, 2, 6, nmax))

    for dlist in dlists:
        (k1,k2,k3,k4) = dlist
        V24.append((dlist,[]))

        for k5 in range(max(-nmax+k1+k2+k3+k4,-nmax),
                min(int(floor((k1+k2+k3+k4)/2)),nmax)+1):

            k6 = k1+k2+k3+k4-k5
            clist = (k5,k6)

            if oscEnergy(clist) <= Emax+ tol:
                V22[-1][1].append(clist)


    V24 = Operator(V24, 4, 2, helper)



    V33 = []

    dlists = set()
    for state in basis.stateList:
        dlists.update(gendlists(state, 3, 6, nmax))

    for dlist in dlists:
        (k1,k2,k3) = dlist
        V33.append((dlist,[]))

        for k4 in range(-nmax, nmax+1):

            for k5 in range(max(-nmax+k1+k2+k3-k4,-nmax),
                min(int(floor((k1+k2+k3-k4)/2)),nmax)+1):

                k6 = k1+k2+k3-k4-k5
                clist = (k4,k5,k6)

                if oscEnergy(clist) <= Emax+ tol \
                    and sorted([abs(min(dlist)),abs(max(dlist))]) <= \
                        sorted([abs(min(clist)),abs(max(clist))]):
                            # XXX This is a generalization of the lexicographical
#sorting condition
                    V33[-1][1].append(clist)


    V33 = Operator(V33, 3, 3, helper)

    return V06, V15, V24, V33

