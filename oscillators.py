from profile_support import *
from operator import mul
from functools import reduce
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
import numpy as np
import me
from operators import *


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
        nd: number of annihilation operators
        nc: number of creation operators
        helper: Helper object
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

# List with same lenght as oscList. It contains the total energy of the corresponding
# oscillators
        self.oscEnergies = []

# Combinatorial and phase-space factors of the oscillators
        self.oscFactors = []

# Overall prefactor
        pref = binom(nc+nd,nc)

        clist_pref = helper.clist_pref
        clist_e = helper.clist_e
        clist_count = helper.clist_count

        for i, (dlist,clists) in enumerate(oscillators):

            dlist_e = oscEnergy(dlist)
            dlist_count = Counter(dlist)
            dlist_pref = pref*factorial(nd)/\
                        reduce(mul, (factorial(c)*sqrt(2*omega(n)*L**2)**c for n,c in dlist_count.items()), 1)

            for clist in clists:
                if clist not in clist_pref:
                    clist_e[clist] = oscEnergy(clist)
                    clist_count[clist] = Counter(clist)
                    clist_pref[clist] = factorial(nc)/reduce(mul, (factorial(c)*sqrt(2*omega(n)*L**2)**c for n,c in clist_count[clist].items()), 1)

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
        wavenumber n
        clist: list of creation momenta
        dlist: list of annihilation momenta
        ccount: dictionary {wn: count} where count is the number of instances of wavenumber wn in clist
        dcount: dictionary {wn: count} where count is the number of instances of wavenumber wn in dlist
        """

        wnlist = set(clist+dlist)
        return scipy.array([[n[0],n[1],ccount.get(n,0),dcount.get(n,0)] for n in wnlist], dtype=scipy.int8)


    def computeMatrixElements(self, basis, i, destbasis, ignKeyErr=False):
        """ Compute the matrix elements by applying all the oscillators in the operator
        to an element in the basis
        basis: set of states on which the operator acts
        i: index of the state in the basis
        destbasis: set of states corresponding to the column indices
        ignKeyErr: this must be set to True if the action of an oscillators on an input state
        can generate a state not in destbasis. This applies in the computation of Vhh.
        Otherwise it should be set to False
        """

        return me.computeME(basis, i, destbasis,
                ignKeyErr, self.nd, self.nc, self.dlistPos, self.oscFactors,
                self.oscList, self.oscEnergies)

    def yieldBasis(self, basis, subidx, EL):
        """ Yields a sequence of representation 2 states, by acting with oscillators
        on a subset of states.
        basis: basis of states which are acted upon
        subidx: subset of indices of basis on which to act
        EL: maximal energy of the generated high-energy states
        """


        for idx in subidx:
            for s in me.yieldBasis(basis,idx,EL,self.helper,self.nd,self.nc,\
                    self.dlistPos, self.oscEnergies,self.oscList):
                yield s


def V4OpsHalf(helper, basis=None):
    """ Generate half of the oscillators of the V4 operator
    helper: Helper function of the destination basis
    basis: starting basis. If None, the states are assumed to be all those
    below the energy cutoff defined in helper.
    """

    dlist = ()
    V40 = [(dlist, helper.genMomenta4sets())]
# Generate a LocOperator instance from the computed set of oscillators
    V40 = LocOperator(V40, 0, 4, helper)

    if basis==None:
        allowedWn = helper.allowedWn
    else:
        allowedWn = basis.helper.allowedWn
    # allowedWnList = list(map(lambda x:np.array(x), sorted(allowedWn)))
    allowedWnList = list(map(lambda x:np.array(x), allowedWn))

    V31 = []
    for k1 in allowedWnList:
# The set of annihilation momenta contains just one momentum
        dlist = (tuple(k1),)
        V31.append((dlist, helper.genMomenta3sets(k1)))
    V31 = LocOperator(V31, 1, 3, helper)


    return V40, V31


def V4Ops22(helper, basis=None):
    # XXX Temporary fix. Do not symmetrize for the moment!
    """ Generate the oscillators of the V4 operator with two annihilation and two creation operators.
    helper: Helper function of the destination basis
    basis: starting basis. If None, the states are assumed to be all those
    below the energy cutoff defined in helper.
    """

    V22 = []

    # TODO Sort the pairs with minEnergy function, so that we can break
# the cycle when we go out of the starting basis?
    # allowedWnPairs = list(helper.genMomentaPairs().values())
    allowedWnPairs = helper.genMomentaPairs()

    if basis==None:
        startallowedWnPairs = allowedWnPairs
    else:
        startallowedWnPairs = basis.helper.genMomentaPairs()


    # elist = [list(map(oscEnergy, kpairlist)) for kpairlist in allowedWnPairs]


    # Cycle over total momentum of annihilation operators
    # for wnIdx in range(len(allowedWnPairs)):

        # kpairlist = allowedWnPairs[wnIdx]

        # for i in range(len(kpairlist)):
            # kpair = kpairlist[i]
            # e12 = elist[wnIdx][i]

            # dlist = kpair
            # V22.append((dlist,[]))

            # for j in range(len(kpairlist)):
                # # XXX Need to perforn any checks ?
                # clist = kpairlist[j]
                # V22[-1][1].append(clist)

    for wn, kpairlist1 in startallowedWnPairs.items():

        kpairlist2 = allowedWnPairs[wn]

        for kpair1 in kpairlist1:

            dlist = kpair1
            V22.append((dlist,[]))

            for kpair2 in kpairlist2:
                # XXX Need to perforn any checks ?
                clist = kpair2
                V22[-1][1].append(clist)

    V22 = LocOperator(V22, 2, 2, helper)
    return (V22,)



def minus(t):
    """ Takes the opposite of a tuple """
    return (-t[0],-t[1])


def V2OpsHalf(helper, basis=None):
    """ Generate half of the oscillators of the V2 operator
    helper: helper function of the destination basis.
    basis: starting basis """

    Emax = helper.Emax

    dlist = ()
    V20 = [(dlist, [])]

    allowedWn12 = helper.allowedWn12
    # Select only in half of phase space
    for k1 in allowedWn12:
        k2 = minus(k1)
        clist = (k1,k2)

        V20[-1][1].append(clist)

    V20 = LocOperator(V20, 0, 2, helper)

    V11 = []
    if basis==None:
        allowedWn = helper.allowedWn
    else:
        allowedWn = basis.helper.allowedWn
    # Select in all phase space
    for k1 in allowedWn:
        dlist = (k1,)
        V11.append((dlist,[]))

        k2 = k1
        clist = (k2,)
        V11[-1][1].append(clist)

    V11 = LocOperator(V11, 1, 1, helper)

    return V20, V11
