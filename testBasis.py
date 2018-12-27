import pytest
import random
from phi4 import *
from nlo import *

def test_genbasis2():
    """ Test the generated basis from the vacuum with momentum cutoff """
    m = 1
    L = 5
    Lambda = 5
    ET = 2
    EL = 3*Lambda
    ELp = EL

    bases = Basis.fromScratch(m, L, ET, Lambda)
    bases2 = Basis.fromScratch(m, L, EL, np.inf)

    subidx = {k: [0] for k in (-1,1)}
    Vlist = None
    V22 = None

    E0 = {1:0, -1:m}

    basesH = genHEBases(bases, subidx, EL, ELp)

    for k in (-1,1):
        basisH = basesH[k]

        maxmom = bases2[k].helper.maxmom

        if k==1:
            idxlist = [i for i,s in enumerate(bases2[k].stateList) if occn(s)==4
                    and maxmom(s)<Lambda+tol]
        else:
            helper = bases2[k].helper
            idxlist = [i for i,s in enumerate(bases2[k].stateList) if
                    maxmom(s)<Lambda+tol and (occn(s)==3 or
                    (occn(s)==5 and helper.torepr2(s)[helper.allowedWn[(0,0)]]>0))]

        assert basisH.size == len(idxlist)

        V = genVHl(bases[k], subidx[k], basisH, L)

        a = Phi4(bases2[k])
        Vlist, V22 = a.computePotential(Vlist, V22)

        V2 = subcolumns(subrows(a.V[4], subidx[k]), idxlist)

        prop = 1/(E0[k]-np.array(basisH.energyList))

        deltaE = np.einsum("ij,j,kj", V2.todense(), prop, V2.todense())[0][0]
        deltaE2 = np.einsum("ij,j,kj", V.todense(), prop, V.todense())[0][0]

        assert abs(deltaE-deltaE2)<tol


def test_genbasis():
    """ Test the generated basis from the vacuum """
    m = 1
    L = 5
    Lambda = 5
    ET = 5
    EL = 3*ET
    ELp = EL

    bases = Basis.fromScratch(m, L, ET, Lambda)
    bases2 = Basis.fromScratch(m, L, EL, Lambda)

    subidx = {k: [0] for k in (-1,1)}
    Vlist = None
    V22 = None

    E0 = {1:0, -1:m}

    basesH = genHEBases(bases, subidx, EL, ELp)

    for k in (-1,1):
        basisH = basesH[k]

        if k==1:
            idxlist = [i for i,s in enumerate(bases2[k].stateList) if occn(s)==4]
        else:
            helper = bases2[k].helper
            idxlist = [i for i,s in enumerate(bases2[k].stateList) if occn(s)==3 or
                    (occn(s)==5 and helper.torepr2(s)[helper.allowedWn[(0,0)]]>0)]

        assert basisH.size == len(idxlist)

        V = genVHl(bases[k], subidx[k], basisH, L)

        a = Phi4(bases2[k])
        Vlist, V22 = a.computePotential(Vlist, V22)

        V2 = subcolumns(subrows(a.V[4], subidx[k]), idxlist)

        prop = 1/(E0[k]-np.array(basisH.energyList))

        deltaE = np.einsum("ij,j,kj", V2.todense(), prop, V2.todense())[0][0]
        deltaE2 = np.einsum("ij,j,kj", V.todense(), prop, V.todense())[0][0]

        assert abs(deltaE-deltaE2)<tol


def test_quartic_spec_Lambda():
    Emax = 15
    L = 5
    Lambda = 4
    g2 = 0
    g4 = 1
    vac = -0.115185001425405
    spece = [1.903760713068681, 3.247950933646434, 4.081462591400829]
    speco = [0.91567578388474,  2.958422505762713, 4.336228677868323,
            5.167300737731443]

    Vlist = None
    V22 = None

    bases = Basis.fromScratch(m=1, L=L, Emax=Emax, Lambda=Lambda)
    eigs = {}

    for k in (-1,1):
        a = Phi4(bases[k])
        Vlist, V22 = a.computePotential(Vlist, V22)

        a.setg(0, g2, g4, ct=False)
        a.setmatrix()

        a.computeEigval(neigs=len(speco))
        eigs[k] = a.eigval

    assert abs(eigs[1][0]-vac) < tol
    np.testing.assert_array_almost_equal(eigs[1][1:]-eigs[1][0], spece)
    np.testing.assert_array_almost_equal(eigs[-1]-eigs[1][0], speco)


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


            a1.setmatrix(Emax=Emax2)
            a2.setmatrix()

            a1.computeEigval()
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

            a.setg(0, g2, g4/(factorial(4)), ct=False)
            a.setmatrix()

            a.computeEigval(neigs=len(so))
            eigs[k] = a.eigval

        assert abs(eigs[1][0]-e0) < tol
        np.testing.assert_array_almost_equal(eigs[1][1:]-eigs[1][0], se)
        np.testing.assert_array_almost_equal(eigs[-1]-eigs[1][0], so)



