# The function computeME is the most computationally intensive, and it has been ported to cython 
# Before running the code, this file has to be compiled via:
# python setup.py build_ext --inplace


# cython: linetrace=False

import scipy
import numpy as np
from math import factorial, floor, sqrt
import statefuncs
from statefuncs import Basis, Helper
import itertools
from itertools import combinations, islice, permutations
from scipy import exp, pi
import bisect
from oscillators import *
cimport cython
from cpython cimport array as array
import array

# XXX Warning usign exponential notation sets this to zero!
cdef double tol = 0.00000001

symFactors = [[0. for _ in range(9)] for _ in range(9)]
for ncomp1 in (1,2,4,8):
    for ncomp2 in (1,2,4,8):
        symFactors[ncomp1][ncomp2] = sqrt(ncomp1/ncomp2)


# XXX Review this function, to take Lambda into account?
def filterDlist(dlist, nd, ntot, helper):
    """ Return True if the list of momenta can be annihilated by any oscillators, so that 
    the new state is in the columns basis 
    dlist: list of potential momenta to annihilate
    nd: number of annihilation operators
    ntot: number of creation+annihilation operators
    helper: Helper object
    """
    if nd==ntot:
        return tuple(sum([np.array(d) for d in dlist])) == (0,0)
    elif nd==ntot-1:
        return tuple(sum([np.array(d) for d in dlist])) in helper.allowedWn
    else:
        return True


def gendlists(state, nd, ntot, helper):
    """ Generates a set of all the sublists of momenta in the state that can be annihilated 
    state: input state in representation 1
    nd: number of annihilation operators (number of modes to annihilate)
    ntot: total number of annihilation and creation operators
    helper: Helper object
    """

    x = itertools.chain.from_iterable(([tuple(n)]*Zn for n,Zn in state))
    dlists = set(tuple(y) for y in combinations(x,nd))

    return (dlist for dlist in dlists if filterDlist(dlist, nd, ntot, helper))


# @cython.binding(True)
# @profile
def computeME(basis, i, statePos, ignKeyErr, nd, nc, dlistPos, oscFactors, oscList, oscEnergies):
        """ Compute the matrix elements by applying all the oscillators in the operator
        to an element in the basis
        basis: set of states on which the operator acts
        i: index of the state in the basis (row index)
        statePos: dictionary {state: idx} state is in representation 2 (in tuple form)
            and idx is its position in the column basis 
        ignKeyErr: this must be set to True if the action of an oscillators on an input state
            can generate a state not in the columns basis (e.g. must be True for computing Vhh in the NLO formalism)
        nd: number of annihilation operators 
        nc: number of creation operators
        dlistPos: dictionary {dlist: idx} for fast access of various objects related to a given list of momenta to annihilate
        oscFactors: data structure containing various numerical factors related to the sets of oscillators
        oscList: data structure with lists of oscillators
        oscEnergies: data structure containing energy shifts due to action of operators
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

        helper = basis.helper

        oscEnergy = helper.oscEnergy

        # List of columns indices of generated basis elements
        col = []
        # List of partial matrix elements
        data = []

        # I define these local variables outside the loops for performance reasons
        e = basis.energyList[i]
        # Number of components in symmetry representation of state
        ncomp = basis.ncomp
        ncompi = basis.ncomp[i]
        state = basis.stateList[i]

        statevec = array.array('b', helper.torepr2(state))
        cstatevec = statevec.data.as_chars
        
        allowedWn = helper.allowedWn
        normFactors = helper.normFactors
        Emax = helper.Emax


        # cycle over all the sets of momenta that can be annihilated
        for dlist in gendlists(state, nd, nd+nc, helper):
            try:
                k = dlistPos[dlist]
            except KeyError as err:
                raise err


# Only select the oscillators such that the sum of the state and oscillator energies
# lies within the bounds of the lookupbasis energies
            imin = bisect.bisect_left(oscEnergies[k], 0-e-tol)
            imax = bisect.bisect_left(oscEnergies[k], Emax-e+tol)
            if imax <= imin:
                continue

            oscFactorsSub = array.array('f', oscFactors[k][imin:imax])
            oscListSub = oscList[k][imin:imax]

            for z in range(len(oscListSub)):

                # XXX This is a bit slow. Maybe need to define oscListSub as a C variable?
                osc = oscListSub[z]

                newstatevec = array.copy(statevec)
                cnewstatevec = newstatevec.data.as_chars

                x = oscFactorsSub[z]

                for ii in range(osc.shape[0]):
                    # Index of momentum in representation 2
                    jj = allowedWn[(osc[ii, 0],osc[ii, 1])]
                    Zc = osc[ii, 2]
                    Zd = osc[ii, 3]
                    cnewstatevec[jj] += Zc-Zd
                    x *= normFactors[Zc, Zd, cstatevec[jj]]


                if ignKeyErr:
                    try:
                        j = statePos[bytes(newstatevec)]
                    except KeyError:
                        continue
                else:
                    j = statePos[bytes(newstatevec)]

                x *= symFactors[ncompi][ncomp[j]]

                data.append(x)
                col.append(j)

        return col, data
