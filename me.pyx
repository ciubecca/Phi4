import gc
from sys import getsizeof as sizeof
import scipy, numpy
from math import factorial, floor, sqrt
import statefuncs
from statefuncs import Basis, Helper
import itertools
from itertools import combinations, islice, permutations
from scipy import exp, pi
from scipy.special import binom
import bisect
from oscillators import *
cimport cython
from cpython cimport array as array
import array

# XXX Warning usign exponential notation sets this to zero!
cdef double tol = 0.00000001

# TODO
# parityFactors = [[1, sqrt(2)],[1/sqrt(2),1]]

# TODO Check and vectorize this 
def symFactors(ncomp1, ncomp2):
    return sqrt(ncomp1/ncomp2)


# XXX Review this function, to take Lambda into account?
def filterDlist(dlist, nd, ntot, helper):
    if nd==ntot:
        return tuple(sum([numpy.array(d) for d in dlist])) == (0,0)
    elif nd==ntot-1:
        return tuple(sum([numpy.array(d) for d in dlist])) in helper.allowedWn
    else:
        return True

# def filterDlist(dlist, nd, ntot, helper):
    # nc = ntot-nd
    # ktot = tuple(sum([numpy.array(d) for d in dlist]))

    # if nc == 0:
        # return ktot == (0,0)
# # XXX Check, and fix to take Lambda into account
    # else:
        # return minEnergy(ktot, nc) < Emax
    # return tuple(sum([numpy.array(d) for d in dlist])) == (0,0)
    # elif nd==ntot-1:
        # return tuple(sum([numpy.array(d) for d in dlist])) in allowedWn
    # else:
        # return True


def gendlists(state, nd, ntot, helper):
    """ Generates a list of all the possible combinations of momenta in the state that
    can be annihilated
    state: input state in representation 1
    nd: number of annihilation operators (number of modes to annihilate)
    ntot: total number of annihilation and creation operators
    allowedWn: all the allowed wave numbers in the basis
    """

    x = itertools.chain.from_iterable(([tuple(n)]*Zn for n,Zn in state))
    dlists = set(tuple(y) for y in combinations(x,nd))

    return (dlist for dlist in dlists if filterDlist(dlist, nd, ntot, helper))

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

        helper = basis.helper

        oscEnergy = helper.oscEnergy

        # List of columns indices of generated basis elements
        col = []
        # List of partial matrix elements
        data = []

        # I define these local variables outside the loops for performance reasons
        e = basis.energyList[i]
        # Number of components in symmetry representation of state
        ncomp1 = basis.ncomp[i]
        state = basis.stateList[i]

        statevec = array.array('b', helper.torepr2(state))
        cstatevec = statevec.data.as_chars
        
        # TODO
        # parityList = lookupbasis.parityList
        allowedWn = helper.allowedWn
        normFactors = helper.normFactors
        Emax = helper.Emax

        # cycle over all the sets of momenta that can be annihilated
        for dlist in gendlists(state, nd, nd+nc, helper):

            debug = False
            state0 = [((0,0),6)]
            state1 = [((0, -1), 1), ((0, 1), 1)]
            state2 = [((-1, 0), 1), ((1, 0), 1)]
            # if state==state0 and nd==2 and nc==2:
                # debug = True
                # print("e state:", e)
                # print("dlist:", dlist)

            try:
                k = dlistPos[dlist]
            except KeyError as err:
                raise err


# Only select the oscillators such that the sum of the state and oscillator energies
# lies within the bounds of the lookupbasis energies
            imin = bisect.bisect_left(oscEnergies[k], 0-e-tol)
            imax = bisect.bisect_left(oscEnergies[k], Emax-e+tol)

            if debug:
                print("oscEnergies", oscEnergies[k])
                print("tol", tol)
                print("Emax-e", Emax-e)
                print("Emax-e+tol", Emax-e+tol)
                print("imin, imax", imin, imax)

            if imax <= imin:
                continue

            if debug:
                print("Check 1")

            oscFactorsSub = array.array('f', oscFactors[k][imin:imax])
            oscListSub = oscList[k][imin:imax]

            for z in range(len(oscListSub)):

                osc = oscListSub[z]

                if debug:
                    # print("osc", [list(x) for x in osc])
                    print("osc:", numpy.asarray(osc))

                newstatevec = array.copy(statevec)
                cnewstatevec = newstatevec.data.as_chars

                x = oscFactorsSub[z]

                for ii in range(osc.shape[0]):
                    # Index of momentum in representation 2
                    jj = allowedWn[tuple(osc[ii, 0:2])]
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
    
                # TODO Check and vectorize
                x *= symFactors(ncomp1, basis.ncomp[j])
                data.append(x)
                col.append(j)

        return col, data


