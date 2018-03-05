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
from itertools import groupby
from scipy import exp, pi
from scipy.special import binom
import bisect
import me
from me import *

tol = 0.000000001

def bose(x):
    """ computes the Bose factor of a product of oscillators  """
    return factorial(len(x))/scipy.prod(
        [factorial(sum(1 for _ in group)) for key, group in groupby(x)]
        )


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

        return scipy.array([[n,cosc.get(n,0),dosc.get(n,0)] for n in wnlist],
                dtype=scipy.int8)

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

        return me.computeME(basis, i, lookupbasis, helper, statePos, Erange,
                ignKeyErr, self.nd, self.nc, self.dlistPos, self.oscFactors,
                self.oscList, self.oscEnergies)


    # Generate high energy Hilbert space Hh from low energy Hilbert space Hl
    # as Hh = V*Hl
    def yieldBasis(self, basis, EL):
        """ Return a set of tuples of representation 2 states, all of which are not
        connected by spatial parity transformations.
        basis: set of "selected" low energy states on which to act
        EL: maximal energy of the generated high-energy states
        """

        nmax = self.helper.nmax

        for i, state in enumerate(basis.stateList):

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


# \phi operator. It connects parity odd and parity even sectors
def V1Ops(basis):
    Emax = basis.Emax
    helper = basis.helper

    dlist = ()
    V10 = [((),[(0,)])]
    V10 = LocOperator(V10, 0, 1, helper)

    V01 = [((0,),[])]
    V01 = LocOperator(V01, 1, 0, helper)

    return (V10,V01)


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


def V2OpsSelectedHalf(basis, Emax, idxList=None):
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

    for nd in (0,1):
        nc = 2-nd

        dlists = gendlistsfromBasis(basis, idxList, nmax, nd, 2)
        oscList = []

        for dlist in dlists:
            clists = [clist for clist in createClistsV2(nmax, dlist, nc) if
                    oscEnergy(clist) <= Emax+tol]
            if nd==1:
                clists = [clist for clist in clists if
                    abs(dlist[0])<=abs(clist[0]) ]

            oscList.append((dlist, clists))

        opsList.append(LocOperator(oscList,nd,nc,helper=helper))

    return opsList



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


def V2OpsSelectedFull(basis, Emax, idxList=None):
    """ Selected set of oscillators between some selected low-energy states
    and states with energy <= Emax
    """

    helper = Helper(basis.helper.m, basis.helper.L, max(Emax,basis.Emax))
    nmax = helper.nmax
    oscEnergy = helper.oscEnergy

    if idxList == None:
        idxList = range(basis.size)

    opsList = []

    for nd in (0,1,2):
        nc = 2-nd

        dlists = gendlistsfromBasis(basis, idxList, nmax, nd, 2)
        oscList = []

        for dlist in dlists:
            clists = [clist for clist in createClistsV2(nmax, dlist, nc) if
                    oscEnergy(clist) <= Emax+tol]
            oscList.append((dlist, clists))

        opsList.append(LocOperator(oscList,nd,nc,helper=helper))

    return opsList


def V4OpsSelectedFull(basis, Emax, idxList=None, ndmax=4):
    """ Selected set of oscillators of the full V4 operator between some selected states
    basis: basis which is acted upon
    Emin: minimum energy of the states to be generated
    Emax: maximal energy of states to be generated
    idxList: subset of indices of the basis which is acted upon
    ndmax: maximum number of annihilation operators
    """

    helper = Helper(basis.helper.m, basis.helper.L, max(Emax,basis.Emax))
    nmax = helper.nmax
    oscEnergy = helper.oscEnergy

    if idxList == None:
        idxList = range(basis.size)

    opsList = []

    for nd in range(ndmax+1):
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


