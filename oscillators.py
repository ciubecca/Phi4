import scipy
from math import factorial, floor, sqrt
from statefuncs import Basis
from collections import Counter
import itertools
from itertools import combinations, islice
from scipy import exp, pi, array
import bisect

tol = 10**(-10)


def bose(x):
    """ computes the Bose factor of a product of oscillators  """
    return factorial(len(x))/scipy.prod(list(map(factorial,Counter(x).values())))


# XXX Check
parityFactors = [[1, sqrt(2)],[1/sqrt(2),1]]


class Operator():
    """
    Collection of oscillators with fixed number of creation and annihilation operators
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

        self.dlistPos = {}
        self.oscList = []
        self.oscEnergies = []
        self.oscFactors = []

        def f(clist, dlist):
            ret = []
            wnlist = set(clist+dlist)
            cosc = Counter(clist)
            dosc = Counter(dlist)
            for n in wnlist:
                ret.append((n, cosc.get(n,0), dosc.get(n,0)))
            return ret


        for i, (dlist,clists) in enumerate(oscillators):
            clists = list(sorted(clists,key=helper.oscEnergy))

            self.dlistPos[dlist] = i

            self.oscList.append([f(clist,dlist) for clist in clists])

            self.oscEnergies.append([helper.oscEnergy(clist)-helper.oscEnergy(dlist)
                for clist in clists])

            self.oscFactors.append([bose(clist)*bose(dlist)*scipy.special.binom(nc+nd,nc)\
                    *scipy.prod([1/sqrt(2*omega(n)*L) for n in clist+dlist])
                    for clist in clists])



    # @profile
    def computeMatrixElements(self, basis, i, lookupbasis):
        # List of columns indices of generated basis elements
        col = []
        # List of partial matrix elements
        data = []

        # I define these local variables outside the loops for performance reasons
        e = basis.energyList[i]
        p = basis.parityList[i]
        state = basis.stateList[i]
        statevec = self.helper.torepr2(state)
        statePos = lookupbasis.statePos
        parityList = lookupbasis.parityList
        Emin = lookupbasis.Emin
        Emax = lookupbasis.Emax
        # print("state", state)
        # print("Emin", Emin)
        # print("Emax", Emax)
        nmax = self.helper.nmax

        normFactors = self.helper.normFactors

        # List of all the groups of momenta of the state that can be annihilated
        for dlist in basis.stateDlists[self.nd][i]:

            k = self.dlistPos[dlist]

            imin = bisect.bisect_left(self.oscEnergies[k], Emin-e-tol)
            imax = bisect.bisect_left(self.oscEnergies[k], Emax-e+tol)

            if imax <= imin:
                continue

            # print("imin, imax", imin, imax)

            # print(self.oscEnergies[k][imin], self.oscEnergies[k][imax-1])

            # imax = bisect.bisect_left(self.oscEnergies[k], Emax-e+tol)
            # if imax==0:
                # continue

            # XXX need to copy the array?
            datapart = (self.oscFactors[k][imin:imax])[:]

            # This cycle could be moved to C
            for i, osc in enumerate(self.oscList[k][imin:imax]):
                newstatevec = statevec[:]

                for n,Zc,Zd in osc:
                    newstatevec[n+nmax] += Zc-Zd
                    datapart[i] *= normFactors[Zc, Zd, statevec[n+nmax]]
                j = statePos[tuple(newstatevec)]

                datapart[i] *= parityFactors[p][parityList[j]]
                col.append(j)

            data += datapart

        return col, data


    # Generate high energy Hilbert space Hh from low energy Hilbert space Hl
    # as Hh = V*Hl
    def genBasis(self, basis, Emin, Emax):

        nmax = self.helper.nmax
        stateset = set()

        for i, state in enumerate(basis):

            statevec = self.helper.torepr2(state)
            e = basis.energyList[i]

            for dlist in basis.stateDlists[self.nd][i]:
                k = self.dlistPos[dlist]

                # XXX Check
                imin = bisect.bisect_left(self.oscEnergies[k], Emin-e+tol)
                imax = bisect.bisect_left(self.oscEnergies[k], Emax-e+tol)
                if imax <= imin:
                    continue

                for i, osc in enumerate(self.oscList[k][imin:imax]):
                    newstatevec = statevec[:]
                    for n,Zc,Zd in osc:
                        newstatevec[n+nmax] += Zc-Zd
                    t1 = tuple(newstatevec)
                    t2 = tuple(newstatevec[::-1])
                    if (t1 not in stateset) and (t2 not in stateset):
                        stateset.add(t1)

        return stateset


def Phi4Operators(helper, basis):

    nmax = max(basis[1].nmax, basis[-1].nmax)
    Emax = max(basis[1].Emax, basis[-1].Emax)

    dlist = ()
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

    V40 = Operator(V40, 0, 4, helper)


    V31 = []
    for k1 in range(-nmax,nmax+1):
        dlist = (k1,)
        V31.append((dlist,[]))

        for k2 in range(-nmax,nmax+1):
            for k3 in range(max(-nmax+k1-k2,k2),
                           min(int(floor((k1-k2)/2)),nmax)+1):

                k4 = k1-k2-k3
                clist = (k2,k3,k4)
                # NOTE The check on dlist is useless here but it'ŝ needed
                # if we generalize the code to other operators
                if helper.oscEnergy(clist) <= Emax+tol\
                    and helper.oscEnergy(dlist) <= Emax+tol:
                    V31[-1][1].append(clist)

    V31 = Operator(V31, 1, 3, helper)


    V22 = []
    for k1 in range(-nmax,nmax+1):
        for k2 in range(k1,nmax+1):
            dlist = (k1,k2)
            V22.append((dlist,[]))

            for k3 in range(max(-nmax+k1+k2,-nmax),
                    min(int(floor((k1+k2)/2)),nmax)+1):

                k4 = k1+k2-k3
                clist = (k3,k4)

                # NOTE The check on dlist is useless here but it'ŝ needed
                # if we generalize the code to other operators
                if helper.oscEnergy(dlist) <= Emax+tol and\
                    helper.oscEnergy(clist) <= Emax+tol and\
                    sorted([abs(k3),abs(k4)])<=sorted([abs(k1),abs(k2)]):
                    # only consider lexicographically ordered part of V22
                    # but also including diagonal part which will be separated below

                    V22[-1][1].append(clist)


    V22 = Operator(V22, 2, 2, helper)

    return V40, V31, V22



def Phi4OperatorsLH(helper, basis, Emin, Emax):

    nmax = helper.nmax

    dlist = ()
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

    V40 = Operator(V40, 0, 4, helper)


    V31 = []
    dlists = set()
    for s in basis.stateDlists[1]:
        dlists.update(s)

    for dlist in dlists:
        k1 = dlist[0]
        V31.append((dlist,[]))

        for k2 in range(-nmax,nmax+1):
            for k3 in range(max(-nmax+k1-k2,k2),
                           min(int(floor((k1-k2)/2)),nmax)+1):

                k4 = k1-k2-k3
                clist = (k2,k3,k4)

                if helper.oscEnergy(clist) <= Emax+tol\
                    and helper.oscEnergy(clist)+tol > helper.oscEnergy(dlist):
                    # We always want to increase the energy
                    V31[-1][1].append(clist)

    V31 = Operator(V31, 1, 3, helper)


    V22 = []
    dlists = set()
    for s in basis.stateDlists[2]:
        dlists.update(s)

    for dlist in dlists:
        (k1,k2) = dlist
        V22.append((dlist,[]))

        for k3 in range(max(-nmax+k1+k2,-nmax),
                min(int(floor((k1+k2)/2)),nmax)+1):

            k4 = k1+k2-k3
            clist = (k3,k4)

            # We always want to increase the energy
            if helper.oscEnergy(clist) >= helper.oscEnergy(dlist) and\
                    helper.oscEnergy(clist) <= Emax+tol:

                V22[-1][1].append(clist)


    V22 = Operator(V22, 2, 2, helper)


    # NOTE V13 and V04 cannot increase the energy

    return V40, V31, V22
