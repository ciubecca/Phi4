import scipy
from profile_support import *
import scipy.sparse.linalg
import scipy.sparse
import math
from math import factorial
from statefuncs import Basis, Helper
from oscillators import *
from operator import attrgetter
import renorm
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

        c = MatrixConstructor(basis)
        # Vlist = {2:V2OpsHalf(basis), 4:V4OpsHalf(basis)}
        Vlist = {2:V2OpsHalf(basis)}
        # for n in (2,4):
        for n in (2,):
            self.V[n] = c.buildMatrix(Vlist[n])*self.L**2
        del c

        self.V[0] = scipy.sparse.eye(basis.size)*self.L**2


    def computeEigval(self, ET, neigs=6):
        """ Compute the eigenvalues for sharp cutoff ET
        neigs: number of eigenvalues to compute
        """

        g0 = self.g0
        g2 = self.g2
        g4 = self.g4

        compH = self.h0 + sum(self.g[n]* self.V[4] for n in (0,2,4))

        # Seed vector
        v0 = scipy.zeros(self.compH.shape[0])
        v0[:10] = 1

        self.eigenvalues = scipy.sort(scipy.sparse.linalg.eigsh(compH, neigs, v0=v0,
                            which='SA', return_eigenvectors=False))

        gc.collect()
