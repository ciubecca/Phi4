import pytest
import random
from phi4 import *

tol = 10**-7

# @pytest.mark.skip(reason="Not needed now")
def test_basis():
    Elist = [10]
    Llist = [10]
    sizelist = [(36179, 36562)]

    for i in range(len(Elist)):
        Emax = Elist[i]
        L = Llist[i]

        bases = Basis.fromScratch(1, L, Emax)

        for j,k in enumerate((1,-1)):
            assert bases[k].size == sizelist[i][j]


def test_sym():
    Elist = [11]
    Llist = [6]

    for i in range(len(Elist)):
        Emax = Elist[i]
        L = Llist[i]

        bases = Basis.fromScratch(1, L, Emax)

        for k in (-1,1):
            a = Phi4(bases[k])
            a.computePotential()
            assert abs(a.V[4]-a.V[4].transpose()).max() < tol

def test_quartic_spec():
    Elist = [6]
    Llist = [5]
    g2 = [0]
    g4 = [24]

    vac = [-0.004076480376270055]
    spece = [[2.035553009079848, 3.199678014214963, 3.356677133229781, 4.082445472190814, 4.21252008745907]]
    speco = [[0.973434499023974, 3.207424273691506, 4.494048612886267, 4.62886587275236, 5.184264427307492, 5.359318131071841]]

    for Emax,L,g2,g4,e0,se,so in zip(Elist,Llist,g2,g4,vac,spece,speco):

        bases = Basis.fromScratch(m=1, L=L, Emax=Emax)
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
