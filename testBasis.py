import pytest
import random
from phi4 import *
from nlo import *

def test_ofpt2():
    """ Test the results of the old fashioned perturbation theory from the
    high energy matrix construction """
    m = 1
    ET = 2
    L = 5
    EL = 17
    Lambda = np.inf

    subidx = {k:[0] for k in (-1,1)}
    E0 = {1:0, -1:m}

    res = {-1:-0.56672212303518199, 1:-0.34285179282771849}

    a = Phi4(m, L, ET, Lambda)
    bases = a.bases
    a.genHEBases(subidx, EL, EL)

    for k in (-1,1):
        basisH = a.basesH[k]
        a.computeHEVs(k)
        prop = 1/(E0[k]-np.array(basisH.energyList))

        idxlist = basisH.subidxlist(EL, Lambda, Emin=2)
        Vsub = subcolumns(a.VHl[k][4], idxlist)
        propsub = prop[idxlist]

        deltaE = np.einsum("ij,j,kj", Vsub.todense(), propsub, Vsub.todense())[0][0]
        np.testing.assert_almost_equal(deltaE, res[k])


def test_genbasis2():
    """ Test the generated basis from the vacuum with momentum cutoff """
    m = 1
    L = 5
    Lambda = 5
    ET = 2
    EL = 3*Lambda
    ELp = EL

    bases = Basis.fromScratch(m, L, ET, Lambda)

    subidx = {k: [0] for k in (-1,1)}

    E0 = {1:0, -1:m}

    basesH = genHEBases(bases, subidx, EL, ELp, V2=False)

    a = Phi4(m, L, ET=EL)
    bases2 = a.bases
    a.computePotential()

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

        V2 = subcolumns(subrows(a.V[k][4], subidx[k]), idxlist)

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

    subidx = {k: [0] for k in (-1,1)}

    E0 = {1:0, -1:m}

    basesH = genHEBases(bases, subidx, EL, ELp, V2=False)

    a = Phi4(m, L, EL, Lambda)
    a.computePotential()
    bases2 = a.bases

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

        V2 = subcolumns(subrows(a.V[k][4], subidx[k]), idxlist)

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

    eigs = {}

    a = Phi4(1, L, Emax, Lambda)
    a.computePotential()

    for k in (-1,1):

        a.setg(0, g2, g4, ct=False)
        a.setmatrix(k)

        a.computeEigval(k, neigs=len(speco))
        eigs[k] = a.eigval[k]

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

        a1 = Phi4(1, L, Emax1, lam)
        a2 = Phi4(1, L, Emax2, lam)

        a1.computePotential()
        a2.computePotential()

        a1.setg(0, g2, g4/(factorial(4)))
        a2.setg(0, g2, g4/(factorial(4)))

        for k in (-1,1):
            a1.setmatrix(k, Emax=Emax2)
            a2.setmatrix(k)

            a1.computeEigval(k)
            a2.computeEigval(k)
            eigs1 = a1.eigval[k]
            eigs2 = a2.eigval[k]

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

        a = Phi4(m, L, Emax, lam)
        a.computePotential()

        for k in (-1,1):
            assert abs(a.V[k][4]-a.V[k][4].transpose()).max() < tol

def test_quartic_spec():
    Elist = [12]
    Llist = [6]
    g2 = [0]
    g4 = [24]

    vac = [-0.1636168421208604]
    spece = [[1.933910936089044, 2.974636035432578, 3.677984551986206,]]
    speco = [[0.933065191544471, 2.992584298393827, 4.050544605726072, 4.715377240194771]]

    for Emax,L,g2,g4,e0,se,so in zip(Elist,Llist,g2,g4,vac,spece,speco):

        # bases = Basis.fromScratch(m=1, L=L, Emax=Emax)
        eigs = {}
        a = Phi4(1, L, Emax)
        a.computePotential()

        for k in (-1,1):

            a.setg(0, g2, g4/(factorial(4)), ct=False)
            a.setmatrix(k)

            a.computeEigval(k, neigs=len(so))
            eigs[k] = a.eigval[k]

        assert abs(eigs[1][0]-e0) < tol
        np.testing.assert_array_almost_equal(eigs[1][1:]-eigs[1][0], se)
        np.testing.assert_array_almost_equal(eigs[-1]-eigs[1][0], so)
