import scipy
from math import factorial, floor, sqrt
from statefuncs import Basis, phi4Info
from collections import Counter
from itertools import combinations, islice
from scipy import exp, pi, array
import bisect

tol = 10**(-10)

def bose(x):
    """ computes the Bose factor of a product of oscillators  """
    return factorial(len(x))/scipy.prod(list(map(factorial,Counter(x).values())))

class Operator():
    """
    Collection of oscillators with fixed number of creation and annihilation operators
    """

    # @profile
    def __init__(self, oscillators, nd, nc, info):
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
        omega = info.omega
        L = info.L
        m = info.m
        self.info = info


        self.dlistPos = {}
        self.doscRepr2 = []
        self.coscRepr2 = []
        self.diffRepr2 = []
        self.oscEnergies = []
        self.oscFactors = []
        self.doscRepr1 = []
        self.coscRepr1 = []
        # self.oscRepr4 = []

        # self.fv = scipy.vectorize(self.info.computeNorm)

        # self.stateDlists = [set(map(lambda x: tuple(sorted(x)),combinations(state,nd)))
                # for state in basis.repr1List]

        for i, (dlist,clists) in enumerate(oscillators):
            clists = list(sorted(clists,key=info.oscEnergy))

            # self.doscRepr1.append(list(Counter(dlist).items()))

            # self.coscRepr1.append([list(Counter(clist).items()) for clist in clists])


            # self.oscRepr4.append()

            self.dlistPos[dlist] = i
            self.doscRepr2.append(info.repr1torepr2(Counter(dlist)))

            self.coscRepr2.append(array([info.repr1torepr2(Counter(clist))
                for clist in clists]))

            self.diffRepr2.append(self.coscRepr2[-1]-self.doscRepr2[-1])

            # XXX Check
            self.oscEnergies.append(array([info.oscEnergy(clist)-info.oscEnergy(dlist)
                for clist in clists]))

            # XXX Check
            self.oscFactors.append(
                    array([bose(clist)*bose(dlist)*scipy.special.binom(nc+nd,nc)\
                    *scipy.prod([1/sqrt(2*omega(n)*L) for n in clist+dlist])
                    for clist in clists]))


    # @profile
    def computeMatrixElements(self, basis, i, lookupbasis):
        # List of columns indices of generated basis elements
        col = []
        # List of partial matrix elements
        data = []

        # I define these local variables outside the loops for performance reasons
        p = basis.parityList[i]
        e = basis.energyList[i]
        state = basis.stateList[i]
        statePos = lookupbasis.statePos
        # print("state", state)
        parityList = lookupbasis.parityList
        Emax = lookupbasis.Emax
        nmax = self.info.nmax
        lookup = lookupbasis.lookup
        # fv = scipy.vectorize(lookup)

        normFactors = self.info.normFactors[:,:,state]
        # print("shape(normFactors)", normFactors.shape)
        # normFactors2 = self.info.normFactors
        parityFactors = self.info.parityFactors[p]

        # List of all the groups of momenta of the state that can be annihilated
        # Call to set() is needed to eliminate duplicates
        # dlists = set(combinations(self.info.repr2torepr1(state), self.nd))

        # print("nd", self.nd)
        # print("dlists", dlists)

        for dlist in basis.stateDlists[self.nd][i]:

            # dlist = tuple(sorted(dlist))
            # print("dlist", dlist)
            k = self.dlistPos[dlist]

            imax = bisect.bisect_left(self.oscEnergies[k], Emax-e+tol)
            if imax==0:
                continue

            oscFactors = self.oscFactors[k][:imax]
            coscRepr2 = array(self.coscRepr2[k][:imax])
            doscRepr2 = array(self.doscRepr2[k])
            diffRepr2 = self.diffRepr2[k][:imax]
            # diffRepr2 = islice(self.diffRepr2[k],imax)

            # coscRepr1 = self.coscRepr1[k][:imax]
            # doscRepr1 = self.doscRepr1[k]

            # print("oscEnergies", self.oscEnergies[k][:imax])
            # print("coscRepr2", coscRepr2)
            # print("doscRepr2", doscRepr2)
            # print("diffRepr2", diffRepr2)

            # coscRepr2 = self.coscRepr2[k][:imax]
            # doscRepr2 = self.doscRepr2[k]

            # x = scipy.full(len(coscRepr1),
                    # scipy.prod([normFactors2[Zn, 0, state[n+nmax]] for n,Zn in doscRepr1]))
            # halfstateList = state - doscRepr2
            # newstateList2 = halfstateList + coscRepr2

            # for j,cosc in enumerate(coscRepr1):
                # x[j] *= scipy.prod([normFactors2[0,Zn,halfstateList[n+nmax]]
                    # for (n,Zn) in cosc])

            # newstateList = state + diffRepr2
            # fv(state+diffRepr2)
            # print("newstateList", newstateList)
            # colpart = [lookup(x) for x in newstateList]
            # colpart = [lookup(state + diff) for diff in diffRepr2]
            colpart = [statePos[tuple(state + diff)] for diff in diffRepr2]
            # print("colpart", colpart)
            # print(lookupbasis[colpart[0]])

            # XXX Can this be optimized?
            # print("normFactors", normFactors)
            # print("shape", normFactors[:,doscRepr2].shape)
            # print("normFactors", normFactors[coscRepr2,doscRepr2])
            # datapart = scipy.prod(normFactors[coscRepr2,doscRepr2], axis=(1,2))


            # print(normFactors[:,doscRepr2].shape)
            x = normFactors[:,doscRepr2].diagonal(axis1=1,axis2=2)
            datapart = x[coscRepr2].diagonal(axis1=1,axis2=2).prod(axis=1)

            # datapart = scipy.prod(self.fv(coscRepr2, doscRepr2, state), axis=1)
            # print("datapart", datapart)
            datapart *= oscFactors*parityFactors[parityList[colpart]]
            # print(datapart)

            col += colpart
            data += datapart.tolist()

        return col, data


    # TODO Generate high energy Hilbert space Hh from low energy Hilbert space Hl
    # as Hh = V*Hl
    def genBasis(self):
        return


# @profile
def Phi4Operators(info):

    Emax = info.Emax
    nmax = info.nmax

    dlist = ()
    V40 = [(dlist, [])]
    for k1 in range(-nmax,nmax+1):
        for k2 in range(k1,nmax+1):
            # NOTE the boundaries for k3 ensure that k3<=k4<=nmax
            for k3 in range(max(-nmax-k1-k2,k2),
                    min(int(floor((-k1-k2)/2)),nmax)+1):

                k4 = -k1-k2-k3
                clist = (k1,k2,k3,k4)

                if info.oscEnergy(clist) <= Emax+tol:
                    V40[-1][1].append(clist)

    V40 = Operator(V40, 0, 4, info)


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
                if info.oscEnergy(clist) <= Emax+tol\
                    and info.oscEnergy(dlist) <= Emax+tol:
                    V31[-1][1].append(clist)

    V31 = Operator(V31, 1, 3, info)


    V22 = []
    for k1 in range(-nmax,nmax+1):
        for k2 in range(k1,nmax+1):
            dlist = (k1,k2)
            V22.append((dlist,[]))

            for k3 in range(max(-nmax+k1+k2,-nmax),
                    min(int(floor((k1+k2)/2)),nmax)+1):

                k4 = k1+k2-k3
                clist = (k3,k4)

                if info.oscEnergy(dlist) <= Emax+tol and\
                    info.oscEnergy(clist) <= Emax+tol and\
                    sorted([abs(k3),abs(k4)])<=sorted([abs(k1),abs(k2)]):
                    # only consider lexicographically ordered part of V22
                    # but also including diagonal part which will be separated below

                    V22[-1][1].append(clist)

                    # print(dlist,clist)

    V22 = Operator(V22, 2, 2, info)

    return V40, V31, V22
