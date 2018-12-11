import pytest
import random
from phi4 import *

tol = 10**-6

# @pytest.mark.skip(reason="Not needed now")
def test_basis():
    Elist = [10]
    Llist = [8]
    sizelist = [(994, 972)]

    for i in range(len(Elist)):
        Emax = Elist[i]
        L = Llist[i]

        bases = Basis.fromScratch(1, L, Emax, sym=True)

        for j,k in enumerate((1,-1)):
            assert bases[k].size == sizelist[i][j]


def test_sym():
    Elist = [11]
    Llist = [6]

    for i in range(len(Elist)):
        Emax = Elist[i]
        L = Llist[i]

        bases = Basis.fromScratch(1, L, Emax, sym=True)

        for k in (-1,1):
            a = Phi4(bases[k])
            a.computePotential()
            assert abs(a.V[4]-a.V[4].transpose()).max() < tol

def test_quartic_spec():
    Elist = [12]
    Llist = [6]
    g2 = [0]
    g4 = [24]

    vac = [-0.1636168421208604]
    spece = [[1.933910936089044, 2.974636035432578, 3.677984551986206,]]
    speco = [[0.933065191544471, 2.992584298393827, 4.050544605726072, 4.715377240194771]]

    for Emax,L,g2,g4,e0,se,so in zip(Elist,Llist,g2,g4,vac,spece,speco):

        bases = Basis.fromScratch(m=1, L=L, Emax=Emax, sym=True)
        eigs = {}

        for k in (-1,1):
            a = Phi4(bases[k])
            a.computePotential()

            a.setg(0, g2, g4/(factorial(4)))

            a.computeEigval(neigs=len(so))
            eigs[k] = a.eigval

        assert abs(eigs[1][0]-e0) < tol
        np.testing.assert_array_almost_equal(eigs[1][1:]-eigs[1][0], se)
        np.testing.assert_array_almost_equal(eigs[-1]-eigs[1][0], so)
