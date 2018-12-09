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
from scipy import exp, pi, array, e, sqrt
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


    def computePotential(self):
        """
        Builds the potential matrices and the free Hamiltonian
        """

        basis = self.basis
        self.h0 = scipy.sparse.spdiags(basis.energyList,0,basis.size,basis.size)

        self.V = {}

        c = MatrixConstructor(basis)
        Vlist = {2:V2OpsHalf(basis), 4:V4OpsHalf(basis)}
        for n in (2,4):
            self.V[n] = c.buildMatrix(Vlist[n])*self.L**2
        del c

        self.V[0] = scipy.sparse.eye(basis.size)*self.L**2

    def setg(self, g0, g2, g4):
        self.g = {}
        self.g[0] = g0
        self.g[2] = g2
        self.g[4] = g4

    def computeEigval(self, neigs=6):
        """ Compute the eigenvalues for sharp cutoff ET
        neigs: number of eigenvalues to compute
        """

        compH = self.h0 + sum([self.g[n]* self.V[n] for n in (0,2,4)])


        # Seed vector
        v0 = scipy.zeros(compH.shape[0])
        v0[:10] = 1

        self.eigval = scipy.sort(scipy.sparse.linalg.eigsh(compH, neigs, v0=v0,
                            which='SA', return_eigenvectors=False))

        gc.collect()
