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
    def __init__(self, L, g, ET):
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

        self.EL = ratioELET*ET
        self.ELp = ratioELpET*ET
        self.ELpp = ratioELppELp*self.ELp
        EL = self.EL
        ELp = self.ELp
        ELpp = self.ELpp

        self.a = phi4.Phi4(self.m, self.L)
        self.a.buildBasis(Emax=self.ET)

        self.a.computePotential(k)

        self.a.setCouplings(g4=self.g)


    def testRaw(self):
        a = self.a

        a.computeEigval(self.k, self.ET, "raw", neigs=self.neigs)
        self.eps = a.eigenvalues["raw"][self.k][0]

        print(self.eps)


    def testLoc(self):
        a = self.a
        k = self.k
        self.a.computeEigval(self.k, self.ET, "renloc", eps=self.eps, neigs=self.neigs)
        self.eps = a.eigenvalues["renloc"][k][0]

        print(self.eps)


    def testTails(self):
        a = self.a
        maxntails = self.maxntails
        k = self.k
        ET = self.ET
        EL = self.EL
        ELp = self.ELp
        ELpp = self.ELpp
        eps = self.eps

        vectorlist = [state for i,state in sorted(enumerate(a.basis[k]), key=lambda x:
                -abs(a.eigenvectors["raw"][k][0][x[0]]))]
        if maxntails != None:
            vectorlist = vectorlist[:maxntails]
        basisl = statefuncs.Basis(k, vectorlist, a.basis[k].helper)

        a.genHEBases(k, basisl, EL=EL, ELpp=ELpp)
        a.computeLEVs(k)
        a.computeHEVs(k)

        a.calcVV3([ELp], eps, test=True)

        a.computeEigval(k, ET, "rentails", EL=EL, ELp=ELp, ELpp=ELpp, eps=eps,
                neigs=self.neigs)
        e = a.eigenvalues["rentails"][k][0]

        print(e)

    def testAll(self):
        self.testRaw()
        self.testLoc()
        self.testTails()

Llist = [5, 8, 10.5]
ETlist = [5, 8, 10.5]
glist = [0.5, 1, 1.5]

for L, ET, g in zip(Llist, ETlist, glist):
    t = testPhi4(L, g, ET)
    t.testAll()
