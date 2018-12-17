import pytest
import random
from phi4 import *

def test_spec_Lambda():
    Elist1 = [14]
    Elist2 = [12]
    Llist = [6]
    lamlist = [4]
    g2 = [0]
    g4 = [24]

    for L,Emax1,Emax2,lam,g2,g4 in zip(Llist,Elist1,Elist2,lamlist,g2,g4):

        # Emax1 > Emax2

        bases1 = Basis.fromScratch(m=1, L=L, Emax=Emax1, Lambda=lam)
        bases2 = Basis.fromScratch(m=1, L=L, Emax=Emax2, Lambda=lam)

        Vlist1 = None
        V221 = None
        Vlist2 = None
        V222 = None

        for k in (-1,1):
            a1 = Phi4(bases1[k])
            a2 = Phi4(bases2[k])

            Vlist1, V221 = a1.computePotential(Vlist1, V221)
            Vlist2, V222 = a2.computePotential(Vlist2, V222)

            a1.setg(0, g2, g4/(factorial(4)))
            a2.setg(0, g2, g4/(factorial(4)))

            a1.computeEigval(Emax=Emax2)
            a2.computeEigval()

            eigs1 = a1.eigval
            eigs2 = a2.eigval

            np.testing.assert_array_almost_equal(eigs1, eigs2)



def test_Lambda():
    Elist = [10,17]
    Llist = [8,5]
    Lamlist = [3, 4]

    m = 1

    for i in range(len(Elist)):
        Emax = Elist[i]
        L = Llist[i]
        lam = Lamlist[i]

        bases1 = Basis.fromScratch(m, L, Emax, lam)
        bases2 = Basis.fromScratch(m, L, Emax)
        bases3 = Basis.fromScratch(m, L, Emax, Emax/2)

        for k in (1,-1):
            maxmom = bases2[k].helper.maxmom
            sl1 = bases1[k].stateList
            sl2 = [s for s in bases2[k].stateList if maxmom(s) <= lam+tol]
            # print(max([maxmom(s) for s in sl2]))

            assert len(sl1) == len(sl2)
            assert len(bases2[k].stateList) == len(bases3[k].stateList)


# @pytest.mark.skip(reason="Not needed now")
def test_occmax():
    """ Test that size of basis built with occmax is equal to subset of full basis """
    Elist = [12]
    Llist = [6]
    occmax = 4
    m = 1

    for i in range(len(Elist)):
        Emax = Elist[i]
        L = Llist[i]

        basesFull = Basis.fromScratch(m, L, Emax, sym=True)
        bases = Basis.fromScratch(m, L, Emax, occmax=occmax, sym=True)

        for j,k in enumerate((1,-1)):
            sub = [s for s in basesFull[k].stateList if basesFull[k].helper.occn(s) <= occmax]
            assert bases[k].size == len(sub)


# @pytest.mark.skip(reason="Not needed now")
def test_basis():
    Elist = [10]
    Llist = [8]
    sizelist = [(994, 972)]

    m = 1

    for i in range(len(Elist)):
        Emax = Elist[i]
        L = Llist[i]

        bases = Basis.fromScratch(m, L, Emax)

        for j,k in enumerate((1,-1)):
            assert bases[k].size == sizelist[i][j]


# @pytest.mark.skip(reason="Not needed now")
def test_sym():
    Elist = [14]
    Llist = [6]
    lamlist = [4]
    m = 1

    for i in range(len(Elist)):
        Emax = Elist[i]
        L = Llist[i]
        lam = lamlist[i]

        bases = Basis.fromScratch(m, L, Emax, lam)

        Vlist = None
        V22 = None

        for k in (-1,1):
            a = Phi4(bases[k])
            Vlist, V22 = a.computePotential(Vlist, V22)
            assert abs(a.V[4]-a.V[4].transpose()).max() < tol

def test_quartic_spec():
    Elist = [12]
    Llist = [6]
    g2 = [0]
    g4 = [24]

    vac = [-0.1636168421208604]
    spece = [[1.933910936089044, 2.974636035432578, 3.677984551986206,]]
    speco = [[0.933065191544471, 2.992584298393827, 4.050544605726072, 4.715377240194771]]

    Vlist = None
    V22 = None

    for Emax,L,g2,g4,e0,se,so in zip(Elist,Llist,g2,g4,vac,spece,speco):

        bases = Basis.fromScratch(m=1, L=L, Emax=Emax)
        eigs = {}

        for k in (-1,1):
            a = Phi4(bases[k])
            Vlist, V22 = a.computePotential(Vlist, V22)

            a.setg(0, g2, g4/(factorial(4)))

            a.computeEigval(neigs=len(so))
            eigs[k] = a.eigval

        assert abs(eigs[1][0]-e0) < tol
        np.testing.assert_array_almost_equal(eigs[1][1:]-eigs[1][0], se)
        np.testing.assert_array_almost_equal(eigs[-1]-eigs[1][0], so)
