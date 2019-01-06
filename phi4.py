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
    def __init__(self, m, L, ET, Lambda=np.inf):

        self.m = m
        self.L = L
        self.ET = ET
        self.Lambda = Lambda
        self.bases = Basis.fromScratch(m, L, ET, Lambda)

        self.Vcomp = {}
        self.V = {}
        self.h0comp = {}
        self.h0 = {}
        self.eigval = {}
        self.eigv = {}

        scipy.set_printoptions(precision=15)


    def computePotential(self):
        """
        Builds the potential matrices and the free Hamiltonian
        If Vlist has been computed for another k, reuse it because it is expensive
        """

        Vlist = None
        V22 = None

        for k in (-1,1):
            basis = self.bases[k]
            helper = basis.helper
            self.h0[k] = scipy.sparse.spdiags(basis.energyList,0,basis.size,basis.size)

            self.V[k] = {}

            if Vlist == None:
                Vlist = {2:V2OpsHalf(helper), 4:V4OpsHalf(helper)}
                V22 = V4Ops22(helper)

            for n in (2,4):
                self.V[k][n] = buildMatrix(basis, Vlist[n])*self.L**2
                # XXX Temporary fix
                if n == 4:
                    self.V[k][n] += buildMatrix(basis, V22,
                            sumTranspose=False)*self.L**2

            self.V[k][0] = scipy.sparse.eye(basis.size)*self.L**2


    def genHEBases(self, tailidx, EL, ELp):

        self.EL = EL
        self.ELp = ELp
        self.tailidx = tailidx
        self.basesH = genHEBases(self.bases, tailidx, EL, ELp)

        self.Vcomp = {}
        self.h0 = {}


    def computeHEVs(self, k):
        """
        Compute the matrices involving the high-energy states below EL
        """

        # Matrix subscript notation:
        # "l": selected low-energy state with energy <= ET
        # "L": generic low-energy state with energy <= ET
        # "h": selected high-energy state with energy <= EL'
        # "H": selected high-energy states with energy <= EL

        basis = self.bases[k]
        subidx = self.tailidx[k]
        basisH = self.basisH[k]
        helperH = basisH.helper

        #################################
        # Generate the VlH matrices
        #################################

        print("Computing VHl...")

        self.VHl[k] = {}

        for (n,Vops) in ((2,V2OpsSelectedFull(basis,helperH)),
                        (4,V4OpsSelectedFull(basis,helperH))):
### TODO pass energy max(ELp, EL) ?
            Vlist = VOps(basis, basisH.helper, subidx)
            self.VHl[k][n] = buildMatrix(basis, Vlist, destbasis=basisH,
                    idxList=subidx, sumTranspose=False)*L**2


        ##############################
        # Generate the VHL matrix
        ##############################

        print("Computing VHL...")

        self.VHL[k] = {}
        self.VLH[k] = {}

### XXX Should I use the full V operator instead of generating it from
### the basis?
        # for (n,Vops) in ((2,V2OpsSelectedFull), (4,V4OpsSelectedFull)):


        for (n,Vops) in ((2,V2OpsHalf), (4,V4OpsHalf), (4,V22OpsHalf)):
        # for (n,Vops) in ((2,V2OpsSelectedFull), (4,V4OpsSelectedFull)):

### TODO pass energy EL
            Vlist = VOps(basisH.helper, basis)
            self.VHL[k][n] = buildMatrix(basis, Vlist, destbasis=basisH,
                    ignKeyErr=True, sumTranspose=False)*L**2
            self.VLH[k][n] = self.VHL[k][n].transpose()


    def setg(self, g0, g2, g4, ct=True):
        """ ct: add logarithmic mass counterterm """
        self.g = {}
        Lambda = self.Lambda
        m = self.m

        if ct:
# The counterterm was computed by defining the Hamiltonian as g2 V2  + g4/(4 !) V4.
# Instead in the code g4 is not divided by 4!
# XXX The o(1) term is not correct
            dg2 = -(g4*factorial(4))**2*1/(12*(4*pi)**2)*log(Lambda/m)
        else:
            dg2 = 0

        self.g[0] = g0
        self.g[2] = g2 - dg2
        self.g[4] = g4



    def setmatrix(self, k, Emax=np.inf, Lambda=np.inf):

        basis = self.bases[k]

        if Emax<basis.Emax-tol or Lambda<basis.Lambda-tol:
            subidx = basis.subidxlist(Emax, Lambda)
            self.Vcomp[k] = {n: submatrix(self.V[k][n], subidx) for n in (0,2,4)}
            self.h0comp[k] = submatrix(self.h0[k], subidx)
        else:
            self.Vcomp[k] = self.V[k]
            self.h0comp[k] = self.h0[k]


    def computeEigval(self, k, neigs=6, eigv=False):
        """ Compute the eigenvalues for sharp cutoff ET
        neigs: number of eigenvalues to compute
        eigv: return eigenvectors
        """

        compH = self.h0comp[k] + sum([self.g[n]*self.Vcomp[k][n] for n in (0,2,4)])

        # Seed vector
        v0 = scipy.zeros(compH.shape[0])
        v0[:10] = 1

        if eigv:
# XXX Check this returns them sorted
            self.eigval[k], self.eigv[k] = scipy.sparse.linalg.eigsh(compH,
                neigs, v0=v0, which='SA', return_eigenvectors=True)
        else:
            self.eigval[k] = scipy.sort(scipy.sparse.linalg.eigsh(compH,
                neigs, v0=v0, which='SA', return_eigenvectors=False))

        gc.collect()
