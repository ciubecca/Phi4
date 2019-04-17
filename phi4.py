# This is one of the main interfaces of the program
# Contains a class to construct the matrices and compute eigenvalues


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
from nlo import *
from counterterms import *
from scipy.integrate import quad
from scipy import exp, pi, array, e, sqrt, log
import numpy as np


class Phi4():
    """ main class """
    def __init__(self, m, L, ET, Lambda=np.inf, momcut=True):
        """
        m: mass
        L: torus side
        ET: energy cutoff
        Lambda: momentum cutoff
        momcut: wether to use the momentum cutoff
        """

        self.momcut = momcut
        self.m = m
        self.L = L
        self.ET = ET
        self.Lambda = Lambda
        self.bases = Basis.fromScratch(m, L, ET, Lambda)

        self.Vcomp = {}
        self.V = {}
        self.h0comp = {}
        self.nonlocct = {}
        self.h0 = {}
        self.eigval = {}
        self.eigv = {}

        scipy.set_printoptions(precision=15)


    def computePotential(self):
        """
        Builds the potential matrices and the free Hamiltonian
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


    # NOTE Not implemented yet
    def genHEBases(self, tailidx, EL, ELp):
        """ Generate a set of high energy states from a set of "tails". It can be used
        both for old fashioned PT computations and for a future implementation of the NLO formalism
        tailidx: set of indices of the tail states
        EL: maximal energy of the high energy basis, for the VLH operator
        ELp: maximal energy of the high energy basis, for the VlH operator
        """

        self.EL = EL
        self.ELp = ELp
        self.tailidx = tailidx
        self.basesH = genHEBases(self.bases, tailidx, EL, ELp)

        self.Vcomp = {}
        self.h0 = {}


    # NOTE Not implemented yet
    def computeHEVs(self, k):
        """
        Compute the matrices involving the high-energy states below EL
        k: parity quantum number
        """

        # Matrix subscript notation:
        # "l": selected low-energy state with energy <= ET
        # "L": generic low-energy state with energy <= ET
        # "h": selected high-energy state with energy <= EL'
        # "H": selected high-energy states with energy <= EL

        self.VHL = {}
        self.VLH = {}
        self.VHl = {}
        L = self.L

        for k in (-1,1):

            basis = self.bases[k]
            subidx = self.tailidx[k]
            basisH = self.basesH[k]
            helperH = basisH.helper

            #################################
            # Generate the VlH matrices
            #################################

            print("Computing VHl...")

            self.VHl[k] = {}

            for (n,VOps) in ((2, V2OpsSelectedFull), (4, V4OpsSelectedFull)):
### TODO pass energy max(ELp, EL) ?
                # Vlist = VOps(basis, basisH.helper, subidx, half=True)
                # XXX Is this efficient? We are not restricting to half the matrix
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

# NOTE we only need half because we don't need the decrease the energy
            for (n,VOps) in ((2,V2OpsHalf), (4,V4OpsHalf), (4,V4Ops22)):

### TODO pass energy EL
                Vlist = VOps(basisH.helper, basis)
                self.VHL[k][n] = buildMatrix(basis, Vlist, destbasis=basisH,
                        ignKeyErr=True, sumTranspose=False)*L**2
                self.VLH[k][n] = self.VHL[k][n].transpose()


    # XXX Not tested
    def computeVhh(self, k):

        print("Computing Vhh")

        basisH = self.basesH[k]
        helperH = basisH.helper

        for k in (1,):
            # Add both V2 and V4 corrections
            for (Vops,sym) in ((4,V4OpsSelectedHalf,True), (4,V4OpsSelected22,True)):
                Vlist = VOps(basisH, helperH)
                self.Vhh[k] = buildMatrix(basisH, Vlist, ignKeyErr=False, sumTranspose=sym)*self.L**2



    def setg(self, g0, g2, g4, ct=True, cutoff=None, impr=False):
        """
        Set the values of the coupling constants
        g0, g2, g4: Coupling constants
        ct: add logarithmic mass counterterm
        cutoff: value of the cutoff (will be either Lambda or ET depending on self.momcut
        impr: add improvement terms
        """

        self.g = {}
        m = self.m

        if impr==True:
            raise ValueError("Not implemented yet")

        dg2 = 0

        if ct:
            if self.momcut:
                if cutoff==None:
                    cutoff = self.Lambda

# The counterterm was computed by defining the Hamiltonian as g2 V2  + g4/(4 !) V4.
# Instead in the code g4 is not divided by 4!
# NOTE Cancel leading perturbative mass correction exactly
                dg2 = -g4**2*ct2Lam(cutoff, m)
                # XXX The o(1) term is not correct
                # dg2 = (g4*factorial(4))**2*1/(12*(4*pi)**2)*log(Lambda/m)

            else:
                if cutoff==None:
                    cutoff = self.ET

                dg2 = -(g4)**2*ct2ET(cutoff, m)

        self.g[0] = g0
        self.g[2] = g2 + dg2
        self.g[4] = g4



    def setmatrix(self, k, Emax=None, Lambda=None):
        """ Set the matrices used in defining the Hamiltonian
        k: parity quantum number
        Emax: actual energy cutoff to be used in the computation. If smaller than Emax defined in this class, submatrices will be taken.
            This option is used to save time when computing eigenvalues at different values of Emax
        Lambda: same as above, for the momentum cutoff
        """

        m = self.m
        basis = self.bases[k]

        if Emax != None or Lambda != None:
            if Emax==None:
                Emax = self.Emax
            if Lambda==None:
                Lambda = self.Lambda
            subidx = basis.subidxlist(Emax, Lambda)
            self.Vcomp[k] = {n: submatrix(self.V[k][n], subidx) for n in (0,2,4)}
            self.h0comp[k] = submatrix(self.h0[k], subidx)

        else:
            self.Vcomp[k] = self.V[k]
            self.h0comp[k] = self.h0[k]

        elist = self.h0comp[k].diagonal()
        if self.momcut == False:
            if Emax==None:
                Emax = self.Emax

            self.nonlocct[k] = -scipy.sparse.diags([[ct0ETnonloc(Emax,e,m) for e in elist]],[0])
        else:
            self.nonlocct[k] = scipy.sparse.diags([[0 for e in elist]],[0])


    def computeEigval(self, k, neigs=6, eigv=False, nonlocct=True):
        """ Compute the eigenvalues for sharp cutoff ET
        k: parity quantum number
        neigs: number of eigenvalues to compute
        eigv: return eigenvectors
        nonlocct: add the non-local counterterm
        """

        compH = self.h0comp[k] + sum([self.g[n]*self.Vcomp[k][n] for n in (0,2,4)])

        if nonlocct:
            compH += self.L**2*self.g[4]**2*self.nonlocct[k]


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
