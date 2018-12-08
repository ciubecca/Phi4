import gc
from sys import getsizeof as sizeof
import scipy, numpy
from math import factorial, floor, sqrt
import statefuncs
from statefuncs import Basis
import itertools
from statefuncs import Helper
from itertools import combinations, islice, permutations
from scipy import exp, pi
from scipy.special import binom
import bisect
from oscillators import *
cimport cython
from cpython cimport array as array
import array

cdef double tol = 10**(-10)

# TODO
# parityFactors = [[1, sqrt(2)],[1/sqrt(2),1]]


def filterDlist(dlist, nd, ntot, allowedWn):
    if nd==ntot:
        # FIXME To perform the sum we need to convert the tuples of momenta to arrays
        return tuple(sum(dlist)) == (0,0)
    elif nd==ntot-1:
        return tuple(sum(dlist)) in allowedWn
    else:
        return True


def gendlists(state, nd, ntot, allowedWn):
    """ Generates a list of all the possible combinations of momenta in the state that
    can be annihilated
    state: input state in representation 1
    nd: number of annihilation operators (number of modes to annihilate)
    ntot: total number of annihilation and creation operators
    allowedWn: all the allowed wave numbers in the basis
    """

    x = itertools.chain.from_iterable(([n]*Zn for n,Zn in state))
    dlists = set(tuple(y) for y in combinations(x,nd))
    return (dlist for dlist in dlists if filterDlist(dlist, nd, ntot, allowedWn))


def computeME(basis, i, statePos, ignKeyErr, nd, nc, dlistPos, oscFactors, oscList, oscEnergies):
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

        cdef double x
        cdef array.array statevec, newstatevec
        cdef char *cstatevec
        cdef char *cnewstatevec
        cdef char[:,:] osc
        cdef float[:] oscFactorsSub
        cdef char Zc, Zd, nmax
        cdef int z, ii, jj
        cdef double[:,:,:] normFactors

        # List of columns indices of generated basis elements
        col = []
        # List of partial matrix elements
        data = []

        # I define these local variables outside the loops for performance reasons
        e = basis.energyList[i]
        # TODO
        # p = basis.parityList[i]
        state = basis.stateList[i]

        statevec = array.array('b', helper.torepr2(state))
        cstatevec = statevec.data.as_chars
        
        # TODO
        # parityList = lookupbasis.parityList
        allowedWn = helper.allowedWn
        normFactors = helper.normFactors
        Emax = helper.Emax

        # cycle over all the sets of momenta that can be annihilated
        for dlist in gendlists(state, nd, nd+nc, allowedWn):

            k = dlistPos[dlist]

# Only select the oscillators such that the sum of the state and oscillator energies
# lies within the bounds of the lookupbasis energies
            imin = bisect.bisect_left(oscEnergies[k], 0-e-tol)
            imax = bisect.bisect_left(oscEnergies[k], Emax-e+tol)

            if imax <= imin:
                continue

            oscFactorsSub = array.array('f', oscFactors[k][imin:imax])
            oscListSub = oscList[k][imin:imax]

            for z in range(len(oscListSub)):
                osc = oscListSub[z]

                newstatevec = array.copy(statevec)
                cnewstatevec = newstatevec.data.as_chars

                x = oscFactorsSub[z]

                for ii in range(osc.shape[0]):
                    Zc = osc[ii, 1]
                    Zd = osc[ii, 2]
                    # Index of momentum in representation 2
                    jj = allowedWn[osc[ii, 0]]
                    cnewstatevec[jj] += Zc-Zd
                    x *= normFactors[Zc, Zd, cstatevec[jj]]

                if ignKeyErr:
                    try:
                        j = statePos[bytes(newstatevec)]
                    except KeyError:
                        continue
                else:
                    j = statePos[bytes(newstatevec)]
    
                # TODO
                # x *= parityFactors[p][parityList[j]]
                data.append(x)
                col.append(j)

        return col, data


