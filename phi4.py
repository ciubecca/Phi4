import scipy
from profile_support import *
import scipy.sparse.linalg
import scipy.sparse
import math
from math import factorial
from statefuncs import Basis, Helper
from oscillators import *
from operator import attrgetter
import gc
from matrix import *
from scipy.integrate import quad
from scipy import exp, pi, array, e, sqrt, log
from sys import getsizeof as sizeof
import numpy as np


class Phi4():
    """ main class """
    def __init__(self, basis):

        helper = basis.helper
        self.m = helper.m
        self.L = helper.L
        self.k = basis.k
        self.basis = basis

        scipy.set_printoptions(precision=15)


    def computePotential(self, Vlist=None, V22=None):
        """
        Builds the potential matrices and the free Hamiltonian
        If Vlist has been computed for another k, reuse it because it is expensive
        """

        basis = self.basis
        helper = basis.helper
        self.h0 = scipy.sparse.spdiags(basis.energyList,0,basis.size,basis.size)

        self.V = {}

        c = MatrixConstructor(basis)

        if Vlist == None:
            Vlist = {2:V2OpsHalf(helper), 4:V4OpsHalf(helper)}
            V22 = V4Ops22(helper)

        for n in (2,4):
            self.V[n] = c.buildMatrix(Vlist[n])*self.L**2
            # XXX Temporary fix
            if n == 4:
                self.V[n] += c.buildMatrix(V22, sumTranspose=False)*self.L**2
        del c

        self.V[0] = scipy.sparse.eye(basis.size)*self.L**2

        return Vlist, V22



    def setg(self, g0, g2, g4, ct=True):
        self.g = {}
        Lambda = self.basis.Lambda
        m = self.m

        if ct:
# The counterterm was computed by defining the Hamiltonian as g2 V2  + g4/(4 !) V4.
# Instead in the code g4 is not divided by 4!
            dg2 = -(g4*factorial(4))**2*1/(12*(4*pi)**2)*log(Lambda/m)
        else:
            dg2 = 0

        self.g[0] = g0
        self.g[2] = g2 - dg2
        self.g[4] = g4

    def setmatrix(self, Emax=np.inf, Lambda=np.inf):

        if Emax<self.basis.Emax-tol or Lambda<self.basis.Lambda-tol:
            subidx = self.basis.subidxlist(Emax, Lambda)
            self.Vcomp = {n: submatrix(self.V[n], subidx) for n in (0,2,4)}
            self.h0comp = submatrix(self.h0, subidx)
        else:
            self.Vcomp = self.V
            self.h0comp = self.h0


    def computeEigval(self, neigs=6):
        """ Compute the eigenvalues for sharp cutoff ET
        neigs: number of eigenvalues to compute
        """

        compH = self.h0comp + sum([self.g[n]*self.Vcomp[n] for n in (0,2,4)])

        # Seed vector
        v0 = scipy.zeros(compH.shape[0])
        v0[:10] = 1

        self.eigval = scipy.sort(scipy.sparse.linalg.eigsh(compH, neigs, v0=v0,
                            which='SA', return_eigenvectors=False))

        gc.collect()
