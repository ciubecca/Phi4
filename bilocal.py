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

tol = 0.000000001

# XXX Check
parityFactors = [[1, sqrt(2)],[1/sqrt(2),1]]

def pickN(state, N, occ):
    """ Pick all possible subset of particles of a given state """
    if N==0:
        return [[]]
    n,Zn = state[0]

    r = range(max(0,N-occ+Zn), min(Zn+1,N+1))
    return [[(n,Zm)]+x for Zm in r for x in pickN(state[1:],N-Zm,occ-Zn)]


def pickNM(state, N, M, occ):
    """ Divide a state into two subsets of particles """
    if state == []:
        return [[]]
    n,Zn = state[0]
    if N+M == Zn:
        return [[(n,N,M)]]

    r = range(max(0,N-occ+Zn), min(N+1,Zn+1))
    return [[(n,Zm,Zn-Zm)]+x for Zm in r for x in pickNM(state[1:],N-Zm,M-Zn+Zm,occ-Zn)]



def gendlistPairs(state, ndPair, ntotPair, nmax):

    ndTot = sum(ndPair)
    occn = statefuncs.occn(state)
    if occn < ndTot:
        return []

    ret = []

    for dlistTot in pickN(state, ndTot, occn):

        for dlist in pickNM(dlistTot, ndPair[0], ndPair[1], ndTot):

            dlist1 = tuple(itertools.chain.from_iterable(([n]*Zn1
                for n,Zn1,Zn2 in dlist)))
            dlist2 = tuple(itertools.chain.from_iterable(([n]*Zn2
                for n,Zn1,Zn2 in dlist)))

            if filterDlist(dlist1,ndPair[0],ntotPair[0],nmax) and\
                filterDlist(dlist2,ndPair[1],ntotPair[1],nmax):
                ret.append((dlist1, dlist2))

    return ret


class BilocOperator():
    def __init__(self, JointOscList, ndPair, ncPair, helper):
        """
        oscillators: list of tuples. The first element of the tuple is a tuple of
        wavenumbers of annihilation operators, and the second element a list of
        tuples of wavenumbers of creation operators
        basis: basis on which the operator will act
        nd: number of annihilation operators
        nc: number of creation operators
        """

        self.ndPair = ndPair
        self.ncPair = ncPair
        self.ntotPair = tuple((nd+nc for nc,nd in zip(ncPair,ndPair)))
        omega = helper.omega
        L = helper.L
        m = helper.m
        self.helper = helper

        self.dlistPairPos = {}
        self.oscList = []
        self.oscEnergies = []
        self.oscFactors = []

        oscEnergy = helper.oscEnergy

        def pairOscEnergy(pair):
            return oscEnergy(pair[0])+oscEnergy(pair[1])

        for i, (dlistPair,clistPairs) in enumerate(JointOscList):
            clistPairs = list(sorted(clistPairs, key=pairOscEnergy))

            self.dlistPairPos[dlistPair] = i

            self.oscList.append([self.torepr1(clistPair,dlistPair)
                for clistPair in clistPairs])

            self.oscEnergies.append([pairOscEnergy(clistPair)-pairOscEnergy(dlistPair)
                for clistPair in clistPairs])

            self.oscFactors.append([
                scipy.prod([bose(clist)*bose(dlist)*binom(nc+nd,nc)*\
                    scipy.prod([1/sqrt(2*omega(n)*L) for n in clist+dlist])
                    for clist,dlist,nc,nd in zip(clistPair,dlistPair,ncPair,ndPair)])
                        for clistPair in clistPairs])


    def torepr1(self, clistPair, dlistPair):
        """ This generates a list of tuples of the form [(n, Zc, Zd),...] from two separate
        tuples of the form (k1,...,kn) and (q1,...,qm), where the k's and q's are respectively
        the creation and annihilation momenta
        Zc and Zd are respectively the number of creation and annihilation operators at
        wavenumber n """
        clist = tuple(sorted(clistPair[0]+clistPair[1]))
        dlist = tuple(sorted(dlistPair[0]+dlistPair[1]))
        wnlist = set(clist+dlist)
        cosc = Counter(clist)
        dosc = Counter(dlist)
        return list((n,cosc.get(n,0),dosc.get(n,0)) for n in wnlist)


    def computeMatrixElements(self, basis, i, lookupbasis, helper, statePos, Erange,
                                ignKeyErr=True):
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
        nmax = helper.nmax
        Emin, Emax = Erange

        normFactors = helper.normFactors

        # cycle over all the sets of momenta that can be annihilated
        for dlistPair in gendlistPairs(state, self.ndPair, self.ntotPair,
                lookupbasis.helper.nmax):

            k = self.dlistPairPos[dlistPair]

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
                        j = statePos[bytes(newstatevec)]
                    except KeyError:
                        continue
                else:
                    j = statePos[bytes(newstatevec)]

                x *= parityFactors[p][parityList[j]]
                data.append(x)
                col.append(j)


        return col, data



def gendlistPairsfromBasis(basis, nmax, ndPair, ntotPair):
    ret = set()

    for state in basis:
        ret.update(gendlistPairs(state=state, ndPair=ndPair,
            ntotPair=ntotPair, nmax=nmax))
    return ret



def V2V4Ops(basis):
    """
    Bilocal operators of :V2 V4:
    """

    helper = basis.helper
    nmax = helper.nmax
    Emax = basis.Emax
    oscEnergy = helper.oscEnergy

    ntotPair = (2,4)

    createClistsV = {2:createClistsV2, 4:createClistsV4}

    opsList = []

    for nd1 in (0,1,2):
        for nd2 in (0,1,2,3,4):
            ndPair = (nd1,nd2)

            ncPair = tuple(ntot-nd for ntot,nd in zip(ntotPair,ndPair))

            dlistPairs = gendlistPairsfromBasis(basis, nmax, ndPair, ntotPair)

            JointOscList = []

            for dlistPair in dlistPairs:

                x1 = createClistsV[2](nmax, dlistPair[0], ncPair[0])
                x2 = createClistsV[4](nmax, dlistPair[1], ncPair[1])


                clistPairs = [(clist1,clist2) for clist1 in x1 for clist2 in x2
                        if oscEnergy(clist1)+oscEnergy(clist2) <= Emax+tol]

                JointOscList.append((dlistPair, clistPairs))

            opsList.append(BilocOperator(JointOscList,ndPair,ncPair,helper=helper))

    return opsList



def V4V4Ops(basis):
    """
    Bilocal operators of :V4 V4:
    """

    helper = basis.helper
    nmax = helper.nmax
    Emax = basis.Emax
    oscEnergy = helper.oscEnergy

    ntotPair = (4,4)

    opsList = []

    for nd1 in (0,1,2,3,4):
        for nd2 in (0,1,2,3,4):
            ndPair = (nd1,nd2)

            ncPair = tuple(ntot-nd for ntot,nd in zip(ntotPair,ndPair))

            dlistPairs = gendlistPairsfromBasis(basis, nmax, ndPair, ntotPair)

            JointOscList = []

            for dlistPair in dlistPairs:

                x1 = createClistsV4(nmax, dlistPair[0], ncPair[0])
                x2 = createClistsV4(nmax, dlistPair[1], ncPair[1])


                clistPairs = [(clist1,clist2) for clist1 in x1 for clist2 in x2
                        if oscEnergy(clist1)+oscEnergy(clist2) <= Emax+tol]

                JointOscList.append((dlistPair, clistPairs))

            opsList.append(BilocOperator(JointOscList,ndPair,ncPair,helper=helper))

    return opsList

