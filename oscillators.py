from profile_support import *
from operator import mul
from functools import reduce
import gc
from sys import getsizeof as sizeof
import scipy
from math import factorial, floor, sqrt
import statefuncs
from statefuncs import Basis, Helper, tol
from collections import Counter
import itertools
from itertools import combinations, islice, permutations
from itertools import groupby
from scipy import exp, pi
from scipy.special import binom
import bisect
import me
import numpy as np
from me import *

debug = False

maxx = 0

# @profile
# XXX Does this have precision issues?
def bose(x):
    """ computes the Bose factor of a product of oscillators  """
    # global maxx
    # if len(x) > maxx:
        # maxx = len(x)
        # print(maxx)
    return factorial(len(x))/scipy.prod(
        [factorial(sum(1 for _ in group)) for key, group in groupby(x)]
        )

# TODO Maybe this is necessary to order the dlists in V4
# def toCanonicalOsclist(osclist):
    # """ Order canonically the momenta in a list of oscillators """
    # return list(sorted(osclist))

class LocOperator():
    """
    Collection of oscillators with fixed number of creation and annihilation operators
    This is convenient to compute matrix elements and generate the high-energy basis
    from a set of tails
    """

    @profile
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
        oscEnergy = helper.oscEnergy
        L = helper.L
        m = helper.m
        self.helper = helper

# This is just a dict where the keys are tuples of annihilation momenta and the
# values are indexes
        self.dlistPos = {}

# This is a list of lists. Each list contains all the possible creation momenta
# corresponding to a given set of annihilation momenta
        self.oscList = []
        self.oscList2 = []

# List with same lenght as oscList. It contains the total energy of the corresponding
# oscillators
        self.oscEnergies = []

# Combinatorial and phase-space factors of the oscillators
        self.oscFactors = []

# Overall prefactor
        pref = binom(nc+nd,nc)

        clist_pref = {}
        clist_e = {}
        clist_count = {}

        for i, (dlist,clists) in enumerate(oscillators):

            dlist_pref = pref*bose(dlist)*reduce(mul,[1/sqrt(2*omega(n)*L**2) for n in dlist],1)
            dlist_e = oscEnergy(dlist)
            dlist_count = Counter(dlist)

            for clist in clists:
                if clist not in clist_pref:
                    clist_pref[clist] = bose(clist)*reduce(mul,[1/sqrt(2*omega(n)*L**2) for n in clist])
                    clist_e[clist] = oscEnergy(clist)
                    clist_count[clist] = Counter(clist)

            clists = list(sorted(clists, key=helper.oscEnergy))

            self.dlistPos[dlist] = i

            # XXX Slow
            self.oscList.append([self.torepr1(clist, dlist, clist_count[clist], dlist_count) for clist in clists])

            self.oscEnergies.append([clist_e[clist]-dlist_e for clist in clists])
            self.oscFactors.append([dlist_pref*clist_pref[clist] for clist in clists])



    def torepr1(self, clist, dlist, ccount, dcount):
        """ This generates a list of tuples of the form [(n, Zc, Zd),...] from two separate
        tuples of the form (k1,...,kn) and (q1,...,qm), where the k's and q's are respectively
        the creation and annihilation momenta
        Zc and Zd are respectively the number of creation and annihilation operators at
        wavenumber n """

        wnlist = set(clist+dlist)
        return scipy.array([[n[0],n[1],ccount.get(n,0),dcount.get(n,0)] for n in wnlist], dtype=scipy.int8)


    def computeMatrixElements(self, basis, i, statePos, ignKeyErr=False):
        """ Compute the matrix elements by applying all the oscillators in the operator
        to an element in the basis
        basis: set of states on which the operator acts
        i: index of the state in the basis
        helper: Helper instance
        statePos: dictionary where the keys are states in representation 2 (in tuple form)
        and the values are their position in the basis
        ignKeyErr: this must be set to True if the action of an oscillators on an input state
        can generate a state not in lookupbasis. This applies only in the computation of Vhh.
        Otherwise it should be set to False
        """

        return me.computeME(basis, i, statePos,
                ignKeyErr, self.nd, self.nc, self.dlistPos, self.oscFactors, self.oscList, self.oscEnergies)


def _genMomentaPairs(helper):
    """ Generate sets of all inequivalent pairs of momenta,
    ordered lexicographically, and indexed by total momentum.
    This is a subroutine used to construct the V22 matrix """

    omega = helper.omega
    minEnergy = helper.minEnergy
    allowedWn = helper.allowedWn
    Emax = helper.Emax

    # Sort 2d momenta lexicographically
    allowedWnList = list(map(lambda x:np.array(x), sorted(allowedWn)))
    l = len(allowedWnList)
    elist = [omega(wn) for wn in allowedWnList]

    allowedWn12 = {}

    for i1 in range(l):
        k1 = allowedWnList[i1]
        e1 = elist[i1]

        for i2 in range(i1,l):
            k2 = allowedWnList[i2]
            k12 = tuple(k1+k2)
            e12 = e1+elist[i2]

            debug = False
            st = tuple(sorted((tuple(k1),tuple(k2))))
            tofind = {((0,-1),(0,1)), ((-1,0), (1,0))}
            # if st in tofind:
                # debug = True
                # print("st", st)

            if k12 not in allowedWn:
                continue

            if debug:
                print("Check1")

            if e12+minEnergy(k12) > Emax+tol:
                continue

            if debug:
                print("Check2")

            if k12 not in allowedWn12:
                allowedWn12[k12] = []

            allowedWn12[k12].append((tuple(k1),tuple(k2)))

    # Sort 2d momenta pairs lexicographically
    return list(map(lambda x: list(sorted(x)), allowedWn12.values()))



@profile
def V4OpsHalf(basis):
    """ Generate half of the oscillators of the V4 operator
    basis: basis of all the low-energy states below ET """

    helper = basis.helper
    omega = helper.omega
    minEnergy = helper.minEnergy
    allowedWn = helper.allowedWn
    Emax = helper.Emax
    oscEnergy = helper.oscEnergy

    # Sort wavenumbers lexicographically, and convert to arrays
    allowedWnList = list(map(lambda x:np.array(x), sorted(allowedWn)))
    l = len(allowedWnList)
    # (wn -> idx) where idx is the position in the ordered list
    allowedWnIdx = {tuple(wn):i for i,wn in enumerate(allowedWnList)}
    elist = [omega(wn) for wn in allowedWnList]

    dlist = ()
# The list of annihilation momenta is empty
    V40 = [(dlist, [])]

    for i1 in range(l):
        k1 = allowedWnList[i1]
        e1 = elist[i1]

        for i2 in range(i1, l):
            k2 = allowedWnList[i2]
            e2 = elist[i2]

            # XXX Check
            if e1+e2+minEnergy(k1+k2, 2) > Emax+tol:
                continue

            for i3 in range(i2,l):
                k3 = allowedWnList[i3]
                e3 = elist[i3]

                k4 = -k1-k2-k3

                if tuple(k4) not in allowedWn:
                    continue

                i4 = allowedWnIdx[tuple(k4)]
                if i4 < i3:
                    continue

                e4 = elist[i4]
                if e1+e2+e3+e4 > Emax+tol:
                    continue

                clist = (tuple(k1),tuple(k2),tuple(k3),tuple(k4))
                V40[-1][1].append(clist)

# Generate a LocOperator instance from the computed set of oscillators
    V40 = LocOperator(V40, 0, 4, helper)



    V31 = []
    for k1 in allowedWnList:
# The set of annihilation momenta contains just one momentum
        dlist = (tuple(k1),)
        # The state must have at least another particle if k1 != 0
        e1 = minEnergy(k1)
        V31.append((dlist,[]))

        for i2 in range(l):
            k2 = allowedWnList[i2]
            e2 = elist[i2]

            # XXX Check
            if e1+e2+minEnergy(k1-k2,2) > Emax+tol:
                continue

            for i3 in range(i2, l):
                k3 = allowedWnList[i3]

                k4 = k1-k2-k3

                if tuple(k4) not in allowedWn:
                    continue

                i4 = allowedWnIdx[tuple(k4)]

                if i4 < i3:
                    continue

                e3 = elist[i3]
                e4 = elist[i4]

                # XXX Check
                if e1+e2+e3+e4 > Emax+tol:
                    continue

                clist = (tuple(k2),tuple(k3),tuple(k4))
                V31[-1][1].append(clist)

    V31 = LocOperator(V31, 1, 3, helper)

    return V40, V31

# @profile
def V4Ops22(basis):
    # XXX Temporary fix
    """ Do not symmetrize for the moment! """

    helper = basis.helper
    omega = helper.omega
    minEnergy = helper.minEnergy
    allowedWn = helper.allowedWn
    Emax = helper.Emax
    oscEnergy = helper.oscEnergy

    V22 = []

    # TODO Change name for this !!
    allowedWn12 = _genMomentaPairs(helper)
    elist = [list(map(oscEnergy, kpairlist)) for kpairlist in allowedWn12]

    # debug = True
    # if debug and basis.k==1:
        # print(allowedWn12[3])

    # Cycle over total momentum of annihilation operators
    # XXX is this a bug?
    for wnIdx in range(len(allowedWn12)):

        kpairlist = allowedWn12[wnIdx]

        debug = False
        if debug and wnIdx == 3 and basis.k==1:
            print(kpairlist)

        for i in range(len(kpairlist)):
            kpair = kpairlist[i]
            e12 = elist[wnIdx][i]

            dlist = kpair
            V22.append((dlist,[]))

            for j in range(len(kpairlist)):
                # XXX Need to perforn any checks ?
                clist = kpairlist[j]
                V22[-1][1].append(clist)

    V22 = LocOperator(V22, 2, 2, helper)
    return (V22,)




# Takes the opposite of a tuple
def minus(t):
    return (-t[0],-t[1])

def V2OpsHalf(basis):
    """ Generate half of the oscillators of the V2 operator
    basis: basis of all the low-energy states below ET """

    helper = basis.helper
    nmax = helper.nmax
    Emax = helper.Emax
    allowedWn12 = helper.allowedWn12
    allowedWn = helper.allowedWn

    dlist = ()
    V20 = [(dlist, [])]

    # Select only in half of phase space
    for k1 in allowedWn12:
        k2 = minus(k1)
        clist = (k1,k2)
        # if debug and basis.k==-1:
            # print("clist:", clist)

        V20[-1][1].append(clist)

    V20 = LocOperator(V20, 0, 2, helper)

    V11 = []
    # Select in all phase space
    for k1 in allowedWn:
        dlist = (k1,)
        V11.append((dlist,[]))

        k2 = k1
        clist = (k2,)
        V11[-1][1].append(clist)

    V11 = LocOperator(V11, 1, 1, helper)

    return V20, V11
