import statefuncs
import phi4
import renorm
import sys
import scipy
import math
import database
import unittest
from numpy import testing as npt

class testPhi4(unittest.TestCase):
    def __init__(self, L, g, ET, eigs):
        self.m = 1
        self.neigs = 1
        self.k = 1
        m = self.m
        k = self.k
        neigs = self.neigs
        self.maxntails = None
        maxntails = self.maxntails
        ratioELET = 3
        ratioELpET = 2
        ratioELppELp = 1.5
        self.L = L
        self.g = g
        self.ET = ET

        self.eigs = eigs

        self.EL = ratioELET*ET
        self.ELp = ratioELpET*ET
        self.ELpp = ratioELppELp*self.ELp
        EL = self.EL
        ELp = self.ELp
        ELpp = self.ELpp

        self.a = phi4.Phi4(self.m, self.L)
        self.a.buildBasis(Emax=self.ET)

        self.a.computePotential(k)

        self.a.setglist(glist=[self.g])


    def testRaw(self):
        a = self.a
        eigs = self.eigs

        a.computeEigval(self.k, self.ET, "raw", neigs=self.neigs)
        self.eps = {self.g: a.eigenvalues[self.g]["raw"][self.k][0]}

        npt.assert_almost_equal(self.eps[self.g],eigs[0])


    def testLoc(self):
        eigs = self.eigs
        a = self.a
        k = self.k
        self.a.computeEigval(self.k, self.ET, "renloc", eps=self.eps, neigs=self.neigs)
        self.eps = {self.g: a.eigenvalues[self.g]["renloc"][k][0]}

        npt.assert_almost_equal(self.eps[self.g], eigs[1])


    def testTails(self):
        eigs = self.eigs
        a = self.a
        k = self.k
        ET = self.ET
        EL = self.EL
        ELp = self.ELp
        ELpp = self.ELpp
        eps = self.eps

        basisl = a.basis[k]

        a.computeLEVs(k, basisl)
        a.genHEBasis(k, EL=EL, ELp=ELp, ELpp=ELpp)
        a.computeHEVs(k)

        a.calcVV3(ELp, eps, test=True)

        a.computeEigval(k, ET, "rentails", EL=EL, ELp=ELp, ELpp=ELpp, eps=eps,
                neigs=self.neigs)
        e = a.eigenvalues[self.g]["rentails"][k][0]

        npt.assert_almost_equal(e,eigs[2], decimal=8)

    def testAll(self):
        self.testRaw()
        self.testLoc()
        self.testTails()

Llist = [5, 8, 10.5, 10, 10]
ETlist = [5, 8, 10.5, 14, 14]
glist = [0.5, 1, 1.5, 2, 2]

eigslist = [
[-0.00268222083952, -0.0609692033023,-0.0117915038649],
[-0.11259005032, -0.348246045516, -0.127496055308],
[-0.437287334374, -0.958064275124, -0.493868803234],
[-0.919579558503, -1.55598186578, -0.922438806401]
]

for L, ET, g, eigs in list(zip(Llist, ETlist, glist, eigslist)):
    print("Testing L={}, ET={}, g={} ...".format(L,ET,g))
    t = testPhi4(L, g, ET, eigs)
    t.testAll()
    print("PASSED")
