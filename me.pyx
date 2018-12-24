import scipy
import numpy as np
from math import factorial, floor, sqrt
import statefuncs
from statefuncs import Basis, Helper
import itertools
from itertools import combinations, islice, permutations
from scipy import exp, pi
import bisect
cimport cython
from cpython cimport array as array
import array
from operators import *

# XXX Warning usign exponential notation sets this to zero!
cdef double tol = 0.00000001

symFactors = [[0. for _ in range(9)] for _ in range(9)]
for ncomp1 in (1,2,4,8):
    for ncomp2 in (1,2,4,8):
        symFactors[ncomp1][ncomp2] = sqrt(ncomp1/ncomp2)



def yieldBasis(basis, idx, EL, helper, nd, nc, dlistPos, oscEnergies, oscList):
    
    stateList = basis.stateList
    torepr2 = helper.torepr2
    allowedWn = helper.allowedWn
    energyList = basis.energyList

    state = stateList[idx]
    e = energyList[idx]
    
    cdef char[:,:] osc
    cdef array.array statevec, newstatevec
    cdef char *cstatevec
    cdef char *cnewstatevec

    statevec = array.array('b', torepr2(state))
    cstatevec = statevec.data.as_chars

    ret = []
    
    for dlist in gendlists(state, nd, nd+nc, helper):

        k =dlistPos[dlist]

        imin = bisect.bisect_left(oscEnergies[k], 0-e-tol)
        imax = bisect.bisect_left(oscEnergies[k], EL-e+tol)

        oscListSub = oscList[k][imin:imax]

        for z in range(len(oscListSub)):

            osc = oscListSub[z]

            newstatevec = array.copy(statevec)
            cnewstatevec = newstatevec.data.as_chars

            for ii in range(osc.shape[0]):
                jj = allowedWn[(osc[ii, 0],osc[ii, 1])]
                Zc = osc[ii, 2]
                Zd = osc[ii, 3]
                cnewstatevec[jj] += Zc-Zd
            
            ret.append(newstatevec)
    return ret



def computeME(basis, i, destbasis, ignKeyErr, nd, nc, dlistPos, 
        oscFactors, oscList, oscEnergies):
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

        helper = destbasis.helper
        statePos = destbasis.statePos

        oscEnergy = helper.oscEnergy

        # List of columns indices of generated basis elements
        col = []
        # List of partial matrix elements
        data = []

        # I define these local variables outside the loops for performance reasons
        e = basis.energyList[i]
        # Number of components in symmetry representation of state
        destncomp = destbasis.ncomp
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

                x *= symFactors[ncompi][destncomp[j]]

                data.append(x)
                col.append(j)

        return col, data


