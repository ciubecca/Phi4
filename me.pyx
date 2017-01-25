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
from scipy import exp, pi
from scipy.special import binom
import bisect
from oscillators import *
cimport cython
from cpython cimport array as carray

cdef double tol = 0.000000001

parityFactors = [[1, sqrt(2)],[1/sqrt(2),1]]

def filterDlist(dlist, nd, ntot, nmax):
    # TODO This can be sped up with the n-SUM algorithm
    if nd==ntot:
        return sum(dlist)==0
    elif nd==ntot-1:
        return abs(sum(dlist))<=nmax
    else:
        return True


# TODO The speed optimization resulting from using the N-Sum algorithm could be important
def gendlists(state, nd, ntot, nmax):
    """ Generates a list of all the possible combinations of momenta in the state that
    can be annihilated
    state: input state in representation 1
    nd: number of annihilation operators (number of modes to annihilate)
    ntot: total number of annihilation and creation operators
    nmax: maximal wavenumber of the "lookup" basis
    """

    x = itertools.chain.from_iterable(([n]*Zn for n,Zn in state))

    dlists = set(tuple(y) for y in combinations(x,nd))
    return (dlist for dlist in dlists if filterDlist(dlist, nd, ntot, nmax))


def computeME(basis, i, lookupbasis, helper, statePos, Erange,
    ignKeyErr, nd, nc, dlistPos, oscFactors, oscList, oscEnergies):
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
        cdef carray.array statevec, newstatevec
        cdef char[:,:] osc
        cdef char n, Zc, Zd
        cdef int ii
        cdef double[:,:,:] normFactors

        # List of columns indices of generated basis elements
        col = []
        # List of partial matrix elements
        data = []

        # I define these local variables outside the loops for performance reasons
        e = basis.energyList[i]
        p = basis.parityList[i]
        state = basis.stateList[i]

        statevec = carray.array('b', helper.torepr2(state))

        parityList = lookupbasis.parityList
        nmax = helper.nmax
        normFactors = helper.normFactors
        Emin, Emax = Erange

        # cycle over all the sets of momenta that can be annihilated
# XXX Check: we replaced lookupbasis.helper.nmax with helper.nmax
        for dlist in gendlists(state, nd, nd+nc, nmax):

            k = dlistPos[dlist]

# Only select the oscillators such that the sum of the state and oscillator energies
# lies within the bounds of the lookupbasis energies
            imin = bisect.bisect_left(oscEnergies[k], Emin-e-tol)
            imax = bisect.bisect_left(oscEnergies[k], Emax-e+tol)

            if imax <= imin:
                continue

            oscFactorsSub = oscFactors[k][imin:imax]
            oscListSub = oscList[k][imin:imax]

            for i in range(len(oscListSub)):
                osc = oscListSub[i]

                newstatevec = carray.copy(statevec)

                x = oscFactorsSub[i]

                for ii in range(osc.shape[0]):
                    n = osc[ii, 0]
                    Zc = osc[ii, 1]
                    Zd = osc[ii, 2]
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


