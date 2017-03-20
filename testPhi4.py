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
        self.decimal = 8
        self.m = 1
        self.neigs = 6
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

        self.a = phi4.Phi4(self.m, self.L, k)
        self.a.buildBasis(Emax=self.ET)

        self.a.computePotential()

        self.a.setglist(glist=[self.g])


    def testRaw(self):
        a = self.a
        eigs = self.eigs

        a.computeEigval(self.ET, "raw", neigs=self.neigs)
        self.eps = {self.g: a.eigenvalues[self.g]["raw"][0]}

        npt.assert_almost_equal(self.eps[self.g],eigs[0], decimal=self.decimal)


    def testLoc(self):
        eigs = self.eigs
        a = self.a
        k = self.k
        self.a.computeEigval(self.ET, "renloc", eps=self.eps, neigs=self.neigs)
        self.eps = {self.g: a.eigenvalues[self.g]["renloc"][0]}

        npt.assert_almost_equal(self.eps[self.g], eigs[1], decimal=self.decimal)


    def testTails(self):
        eigs = self.eigs
        a = self.a
        k = self.k
        ET = self.ET
        EL = self.EL
        ELp = self.ELp
        ELpp = self.ELpp
        eps = self.eps

        basisl = a.basis

        a.computeLEVs(basisl)
        a.genHEBasis(EL=EL, ELp=ELp, ELpp=ELpp)
        a.computeHEVs()

        a.calcVV3(ELp, eps, test=True)

        a.computeEigval(ET, "rentails", EL=EL, ELp=ELp, ELpp=ELpp, eps=eps,
                neigs=self.neigs)
        e = a.eigenvalues[self.g]["rentails"][0]

        npt.assert_almost_equal(e,eigs[2], decimal=self.decimal)

    def testAll(self):
        self.testRaw()
        self.testLoc()
        self.testTails()

Llist = [7, 10]
ETlist = [10, 10]
glist = [1.5, 2.5]

# XXX Change these to values with finite volume corrections
eigslist = [
    [-0.28132186085916155,-0.65692871641605322,-0.29930595576064467],
    [-0.96072791675498348, -2.6579143236349014,-1.0564035979501218]
]

for L, ET, g, eigs in list(zip(Llist, ETlist, glist, eigslist)):
    t = testPhi4(L, g, ET, eigs)
    t.testAll()
