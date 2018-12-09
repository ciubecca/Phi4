from profile_support import *
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
import numpy as np
from me import *

tol = 10**(-10)

# @profile
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

    # @profile
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

# List with same lenght as oscList. It contains the total energy of the corresponding
# oscillators
        self.oscEnergies = []

# Combinatorial and phase-space factors of the oscillators
        self.oscFactors = []

        for i, (dlist,clists) in enumerate(oscillators):
            clists = list(sorted(clists, key=helper.oscEnergy))

            self.dlistPos[dlist] = i

            self.oscList.append([self.torepr1(clist,dlist) for clist in clists])

            self.oscEnergies.append([oscEnergy(clist)-oscEnergy(dlist) for clist in clists])

            self.oscFactors.append([bose(clist)*bose(dlist)*binom(nc+nd,nc)\
                    *scipy.prod([1/sqrt(2*omega(n)*L**2) for n in clist+dlist])
                    for clist in clists])

    # @profile
    def torepr1(self, clist, dlist):
        """ This generates a list of tuples of the form [(n, Zc, Zd),...] from two separate
        tuples of the form (k1,...,kn) and (q1,...,qm), where the k's and q's are respectively
        the creation and annihilation momenta
        Zc and Zd are respectively the number of creation and annihilation operators at
        wavenumber n """

        wnlist = set(clist+dlist)

        cosc = Counter(clist)
        dosc = Counter(dlist)

        # return [[n,cosc.get(n,0),dosc.get(n,0)] for n in wnlist]
        return scipy.array([[n[0],n[1],cosc.get(n,0),dosc.get(n,0)] for n in wnlist], dtype=scipy.int8)


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

# @profile
def V4OpsHalf(basis):
    """ Generate half of the oscillators of the V4 operator
    basis: basis of all the low-energy states below ET """

    helper = basis.helper
    omega = helper.omega
    minEnergy = helper.minEnergy
    allowedWn = helper.allowedWn
    Emax = helper.Emax
    oscEnergy = helper.oscEnergy

    allowedWnList = list(sorted([np.array(wn) for wn in allowedWn], key=omega))
    elist = [omega(wn) for wn in allowedWnList]

    dlist = ()
# The list of annihilation momenta is empty
    V40 = [(dlist, [])]

    for i1, k1 in enumerate(allowedWnList):
        e1 = elist[i1]

        for i2, k2 in enumerate(allowedWnList):
            e2 = elist[i2]

            if e1+e2 > Emax:
                break

            if tuple(k1+k2) not in allowedWn:
                continue

            for i3,k3 in enumerate(allowedWnList):
                e3 = elist[i3]
                etot = e1+e2+e3

                if etot > Emax:
                    break

                ktot = k1+k2+k3
                if tuple(ktot) not in allowedWn or etot+minEnergy(ktot)>Emax:
                    continue

                k4 = -ktot
                clist = tuple(sorted((tuple(k) for k in (k1,k2,k3,k4))))

                if oscEnergy(clist) <= Emax+tol:
                    V40[-1][1].append(clist)

# Generate an LocOperator instance from the computed set of oscillators
    V40 = LocOperator(V40, 0, 4, helper)

    return (V40, )


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
