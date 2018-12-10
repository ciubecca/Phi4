import pytest
import random
from phi4 import *

tol = 10**(-10)

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
    Elist = [12]
    Llist = [6]

    for i in range(len(Elist)):
        Emax = Elist[i]
        L = Llist[i]

        bases = Basis.fromScratch(1, L, Emax)

        for k in (-1,1):
            a = Phi4(bases[k])
            a.computePotential()
            assert abs(a.V[4]-a.V[4].transpose()).max() < tol


