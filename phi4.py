import scipy
import scipy.sparse.linalg
import scipy.sparse
import math
from math import factorial
from statefuncs import Basis, Helper
from oscillators import *
from bilocal import *
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
    def __init__(self, m, L, k):

        self.m = m
        self.L = L
        self.k = k
# Maximum dimension of the chunks for computing Vhh
        self.chunklen = 20000

# VEV of phi^2self.vev = {}

        self.basis = None
        self.h0 = None
        self.V = {}
        self.DeltaH = None
        self.VLH = {}
        self.VHL = {}
        self.VHl = {}
        self.Vll = {}
        self.VlL = {}
        self.VLl = {}
        self.V0V4 = None
        self.V2V4 = None
        self.V4V4 = None
        self.basisH = None
        self.basisl = None
        self.DH3ll = None

        scipy.set_printoptions(precision=15)


    def buildBasis(self, Emax, occmax=None):
        """ Builds the full Hilbert space basis up to cutoff Emax """
        self.basis = Basis.fromScratch(m=self.m, L=self.L, k=self.k, Emax=Emax, occmax=occmax)

    def computePotential(self, other=None):
        """
        Builds the potential matrices and the free Hamiltonian. In the low-energy sector
        other: basis of the opposite basis used to compute the V = \phi matrix element
        """

        basis = self.basis

        self.h0 = scipy.sparse.spdiags(basis.energyList,0,basis.size,basis.size)

        c = MatrixConstructor(basis, basis)
        Vlist = {2:V2OpsHalf(basis), 4:V4OpsHalf(basis)}
        for n in (2,4):
            self.V[n] = c.buildMatrix(Vlist[n], sumTranspose=True)*self.L
        del c

        if other!=None:
            c = MatrixConstructor(self.basis, other)
            Vlist = V1Ops(self.basis)
            self.V[1] = c.buildMatrix(Vlist, sumTranspose=False)*self.L

        # Construct the identity potential matrix
        self.V[0] = scipy.sparse.eye(basis.size)*self.L


    def genHEBasis(self, EL, ELp, ELpp):
        """ Generate a high-energy basis from a set of tails
        k: parity quantum number
        basisl: Basis instance containing the set of tails
        EL: maximal energy of the generated basis for DH2
        ELpp: maximal energy of the generated basis for DH3
        """

        basisl = self.basisl

        self.EL = EL
        self.ELp = ELp
        self.ELpp = ELpp
# "or" in case arguments are None
        Emax = max(EL or 0, ELp or 0, ELpp or 0)

        # Generate all the operators between the selected states and the states
        # in the range [0, Emax]
        # XXX Assuming that V2 does not generate different states
        Vlist = V4OpsSelectedFull(basisl, Emax)
        vectorset = set()

        for V in Vlist:
            for v in V.yieldBasis(basisl, Emax):
                # Don't add twice states connected by parity inversion
                if v not in vectorset and v[::-1] not in vectorset:
                    vectorset.add(v)

        helper = Vlist[0].helper

# Basis of selected states with energy <= Emax. We only need to save
# states in the type 1 representation (most memory consuming) for states
# with energy <= ELp, or ET when ELp=None
        self.basisH = Basis(self.k, vectorset, helper, repr1=False,
                repr1Emax=max(ELp or 0, self.basis.Emax))


# XXX We could compute either Vhl or VHl to save time
    def computeHEVs(self):
        """
        Compute the matrices involving the high-energy states below EL
        """

        # NOTE matrix subscript notation:
        # "l": selected low-energy state with energy <= ET
        # "L": generic low-energy state with energy <= ET
        # "h": selected high-energy state with energy <= EL'
        # "H": selected high-energy states with energy <= EL

        #################################
        # Generate the VlH matrices
        #################################

        print("Computing VHl...")

        basis = self.basisl
        lookupbasis = self.basisH
        Emax = lookupbasis.Emax

        c = MatrixConstructor(basis, lookupbasis)

        Vlist = V4OpsSelectedFull(basis, Emax)
        self.VHl[4] = c.buildMatrix(Vlist)*self.L

        Vlist = V2OpsSelectedFull(basis, Emax)
        self.VHl[2] = c.buildMatrix(Vlist)*self.L

        del c

        ##############################
        # Generate the VHL matrix
        ##############################

        print("Computing VHL...")

        basis = self.basis
        lookupbasis = self.basisH

        # We only need this matrix for DH2, not for DH3
        Emax = self.EL
        Erange = (0,Emax)

        c = MatrixConstructor(basis, lookupbasis, Erange=Erange)

        # We just need the full operator
        Vlist = V4OpsSelectedFull(basis, Emax)
        self.VHL[4] = c.buildMatrix(Vlist, ignKeyErr=True)*self.L
        self.VLH[4] = self.VHL[4].transpose()


        # XXX Assuming that we are adding all the tails
        self.VHL[2] = self.VHl[2]
        self.VLH[2] = self.VHl[2].transpose()

        del c


    def computeLEVs(self):

        ###################################
        # Generate all the "local" matrices on the selected
        # low-energy states
        ###################################

        self.basisl = self.basis

        subOp = SubmatrixOperator.fromSubbasis(self.basis, self.basisl)

        self.Vll = {}
        for n in (0,2,4):
            self.VlL[n] = subOp.subcolumns(self.V[n])
            self.VLl[n] = self.VlL[n].transpose()
            self.Vll[n] = subOp.subrows(self.VlL[n])


        basis = self.basisl
        c = MatrixConstructor(basis, basis)


    def computeDeltaH(self, ren, ET, eps, loc2=True, loc3=True, loc3mix=True,
            nonloc3mix=True, EL=None, ELp=None, ELpp=None, subbasisl=None,
            memdbg=False):
        """
        Compute the full DeltaH = DH2 * (DH2-DH3)^-1 * DH2  matrix
        subbasisl: explicit set of tails used in the computation. If not specified,
        the tails in self.basisl below ET are used
        """

        glist = self.glist

        # These operators allow to take submatrices according to a given energy range
        subOPL = SubmatrixOperator.fromErange(self.basis, (0, ET))

        if ren=="rentails":
            if subbasisl==None:
                subOPl = SubmatrixOperator.fromErange(self.basisl, (0,ET))
            else:
                subOPl = SubmatrixOperator.fromSubbasis(self.basisl, subbasisl)

            self.ntails = len(subOPl.idxList)

            DH2ll, DH2Ll = self.computeDH2(ET=ET, EL=EL, eps=eps, loc2=loc2)

            # Extract the submatrix in the right energy range
            DH2Ll = {g: subOPL.subcolumns(subOPl.subrows(DH2Ll[g])) for g in glist}
            DH2lL = {g: DH2Ll[g].transpose() for g in glist}

            DH2ll = {g: subOPl.sub(DH2ll[g]) for g in glist}

            DH3ll = self.computeDH3(ET=ET, ELp=ELp, ELpp=ELpp, eps=eps,
                        loc3=loc3, loc3mix=loc3mix, nonloc3mix=nonloc3mix,
                        memdbg=memdbg)
            DH3ll = {g: subOPl.sub(DH3ll[g]) for g in glist}

            gc.collect()

            return DH2lL, DH2ll, DH3ll, DH2Ll


        elif ren=="renloc":
            ret = {}
            VLL = {}
            for n in (0,2,4):
# Assume no submatrix is taken
                VLL[n] = self.V[n]

            for g in glist:
                # Dictionary of local renormalization coefficients for the g^2 term

                VV2 = renorm.renVV2(g4=g, g2=self.g[2][g], EL=ET, eps=eps[g]).VV2

                ret[g] = VV2[0]*VLL[0] + VV2[2]*VLL[2] + VV2[4]*VLL[4]

            return ret


    def computeDH2(self, ET, EL, eps, loc2):

        glist=self.glist

        propagatorH = {g: self.basisH.propagator(eps[g], ET, EL) for g in glist}

        #######################################
        # Construct DH2
        #######################################
        # Dictionary of local renormalization coefficients for the g^2 term

        DH2Ll = {}
        DH2ll = {}

        for g in glist:
            g2 = self.g[2][g]
            g4 = g

            VV2 = renorm.renVV2(g4=g, g2=g2, EL=EL, eps=eps[g]).VV2

            MHl = self.VHl[2]*g2 + self.VHl[4]*g4
            MlH = MHl.transpose()
            MLH = self.VLH[2]*g2 + self.VLH[4]*g4

            DH2Ll[g] = MHl*propagatorH[g]*MLH
            DH2ll[g] = MHl*propagatorH[g]*MlH

            if loc2:
                DH2Ll[g] += VV2[0]*self.VLl[0]+VV2[2]*self.VLl[2]+VV2[4]*self.VLl[4]
                DH2ll[g] += VV2[0]*self.Vll[0]+VV2[2]*self.Vll[2]+VV2[4]*self.Vll[4]

        return DH2ll, DH2Ll

    def computeDH3(self, ET, ELp, ELpp, eps, loc3, loc3mix, nonloc3mix,
            memdbg=False):

        glist = self.glist

        print("Computing DH3")

        sizel = self.basisl.size
        DH3ll = {g: scipy.sparse.csc_matrix((sizel, sizel)) for g in glist}

        if memdbg:
            print("memory before computing Vhh", memory_usage())

        VHl = {}
        VlH = {}
        for n in (2,4):
            VHl[n] = self.VHl[n]
            VlH[n] = VHl[n].transpose().tocsc()

        MHl = {}
        MlH = {}
        for g in glist:
            MHl[g] = (self.g[2][g]*VHl[2] + g*VHl[4])
            MlH[g] = MHl[g].transpose().tocsc()


        basis = self.basisH
        # List of basis elements on which we will cycle
        fullIdxList = basis.irange((ET,ELp))
        idxlen = len(fullIdxList)
        chunklen = min(self.chunklen, idxlen)
        idxLists = [fullIdxList[x:x+chunklen] for x in range(0, idxlen, chunklen)]

        print("Total number of chunks: ", len(idxLists))

#########################################################################
# Add the "symmetric" contributions to DH3 by integrating out states with
# energy ET < E < ELp
# In this part we also add the exponentially suppressed term V2
#########################################################################

        # Propagator and projector over the states between ET and ELp
        propagatorh = {g: basis.propagator(eps[g], ET, ELp) for g in glist}

        c = MatrixConstructor(basis, basis, (ET, ELp))

        # Add both V2 and V4 corrections
        for n in (2,4):

            if n==4:
                Vlist = V4OpsSelectedHalf(basis, Emax=ELp, idxList=fullIdxList)
            elif n==2:
                Vlist = V2OpsSelectedHalf(basis, Emax=ELp, idxList=fullIdxList)

            ##############################
            # Generate the Vhh matrix
            ##############################
            for i,idxList in enumerate(idxLists):
                print("Doing chunk", i, "for Vhh")

                VhhHalfPart =  c.buildMatrix(Vlist, ignKeyErr=True,
                        sumTranspose=False, idxList=idxList)*self.L

                VhhDiagPart = scipy.sparse.spdiags(VhhHalfPart.diagonal(),0,
                        basis.size,basis.size).tocsc()

                for g in glist:
                    DH3llPart = MHl[g]*propagatorh[g]*VhhHalfPart*propagatorh[g]\
                            *MlH[g]*self.g[n][g]
                    DH3llPart += DH3llPart.transpose()
                    DH3llPart -= MHl[g]*propagatorh[g]*VhhDiagPart*propagatorh[g]\
                            *MlH[g]*self.g[n][g]
                    DH3ll[g] += DH3llPart

                del VhhHalfPart

            del Vlist
        del c

#########################################################################
# Add the "mixed" contributions to DH3 by integrating out states with
# energy ET < E < ELp on one side and ELp < E < EL on the other
#########################################################################
        if memdbg:
            print("memory before computing VhH", memory_usage())

        if nonloc3mix:

            c = MatrixConstructor(basis, basis, (ELp, ELpp))

            # Propagator and projector over the states between ELp and ELpp
            propagatorH = {g: basis.propagator(eps[g], ELp, ELpp) for g in glist}

            # Add both V2 and V4 corrections
            for n in (2,4):

                if n==4:
                    VHhlist = V4OpsSelectedFull(basis, ELpp, idxList=fullIdxList)
                elif n==2:
                    VHhlist = V2OpsSelectedFull(basis, ELpp, idxList=fullIdxList)

                for i, idxList in enumerate(idxLists):
                    print("doing chunk", i, "for VhH")

                    VHhPart = c.buildMatrix(VHhlist, ignKeyErr=True,
                            sumTranspose=False,idxList=idxList)*self.L

                    for g in glist:
                        DH3llPart = MHl[g]*propagatorh[g]*VHhPart*propagatorH[g]\
                                *MlH[g]*self.g[n][g]
                        DH3llPart += DH3llPart.transpose()
                        DH3ll[g] += DH3llPart

                    del VHhPart

                del VHhlist
            del c

#######################################################################
# Add the "mixed" local-nonlocal contributions to DH3 where some states
# are integrated from ET to ELp while others from ELpp to Infinity
########################################################################

        if loc3mix:
            for g in glist:
                VV2 = renorm.renVV2(g4=g, g2=self.g[2][g], EL=ELpp, eps=eps[g]).VV2

                DH3llPart = VHl[2]*VV2[2]*propagatorh[g]*MlH[g]
                DH3llPart += VHl[4]*VV2[4]*propagatorh[g]*MlH[g]
                DH3llPart += DH3llPart.transpose()
                DH3ll[g] += DH3llPart

#####################################################
# Add the "symmetric local" parts to DH3
#####################################################

        if loc3:
            Vll = {}
            for n in (0,2,4,6):
                Vll[n] = self.Vll[n]

            V0V4 = self.V0V4
            V2V4 = self.V2V4
            V4V4 = self.V4V4

            for g in glist:
                DH3ll[g] += V0V4*self.VV3.V0V4[g]*g**3
                DH3ll[g] += V2V4*self.VV3.V2V4[g]*g**3
                DH3ll[g] += V4V4*self.VV3.V4V4[g]*g**3

                DH3ll[g] += Vll[0]*self.VV3.VV3loc[0][g]*g**3
                DH3ll[g] += Vll[2]*self.VV3.VV3loc[2][g]*g**3
                DH3ll[g] += Vll[4]*self.VV3.VV3loc[4][g]*g**3
                DH3ll[g] += Vll[6]*self.VV3.VV3loc[6][g]*g**3

        return DH3ll


    def calcVV3(self, ELp, eps, test=False):
        print("Calculating VVV renorm coefficients")
        if test == True:
            self.VV3 =  renorm.renVV3Test(glist=self.glist)
        else:
            self.VV3 =  renorm.renVV3(m=self.m, ET=ELp, eps=eps, glist=self.glist)

    def setglist(self, glist):
        L = self.L
        self.glist = glist

        def f(x, L):
            return 1/sqrt(L**2+x**2)*1/(e**sqrt(L**2+x**2)-1)

        def E0(L):
            return -1/(pi*L)*quad(lambda x: x**2*f(x,L), 0, np.inf)[0]

        def z(L):
            return 1/pi*quad(lambda x: f(x,L), 0, np.inf)[0]

        def g2(g, L):
            return 6*g*z(L)

        def g0(g, L):
            return E0(L)/L + 3*z(L)**2*g

        self.g = {}
        self.g[0] = {g : g0(g,L) for g in glist}
        self.g[2] = {g : g2(g,L) for g in glist}
        self.g[4] = {g : g for g in glist}

        self.vev = {g:{} for g in glist}
        self.eigenvalues = {g: {"raw":None, "renloc":None, "rentails":None}
                for g in glist}
        self.eigenvectors = {g: {"raw":None, "renloc":None, "rentails":None}
                for g in glist}


    def computeEigval(self, ET, ren, EL=None, ELp=None, ELpp=None, loc2=True,
            eps=None, neigs=6, subbasisl=None, loc3=True, loc3mix=True,
            nonloc3mix=True, memdbg=False, memsave=True):
        """ Compute the eigenvalues for sharp cutoff ET and local cutoff EL
        ET: ET
        EL: EL. Should be set to None for "raw" and to EL for "renlocal"
        ren: type of eigenvalue. Can be "raw", "renlocal" or "rentails"
        eps: epsilon parameter in the propagator
        neigs: number of eigenvalues to compute
        """

        glist = self.glist

        if ren=="renloc":
            DeltaH = self.computeDeltaH(ET=ET, ren=ren, eps=eps, loc2=loc2)
            compH = {g: self.h0 + g*self.V[4] + self.g[2][g]*self.V[2]
                + self.g[0][g]*self.V[0]  + DeltaH[g] for g in glist}
            self.compSize = compH[glist[0]].shape[0]

            # Seed vector
            v0 = scipy.zeros(self.compSize)
            for i in range(min(10,len(v0))):
                v0[i] = 1.

            for g in glist:
                self.eigenvalues[g][ren], eigenvectorstranspose = \
                    scipy.sparse.linalg.eigsh(compH[g], neigs, v0=v0,
                            which='SA', return_eigenvectors=True)
                    # scipy.sort(scipy.sparse.linalg.eigsh(compH[g], neigs, v0=v0,
                            # which='SA', return_eigenvectors=False))

                self.eigenvectors[g][ren] = eigenvectorstranspose.T

                v = self.eigenvectors[g][ren][0]

                if self.k == 1:
                    v = self.eigenvectors[g][ren][0]
                    v2 = (self.V[2]/self.L).dot(v)
                    norm = np.inner(v,v)
                    np.testing.assert_almost_equal(1,norm)
                    s = np.inner(v2,v)
                    self.vev[g][ren] = s

            return


        if ren=="raw":
            compH = {g: self.h0 + g*self.V[4] + self.g[2][g]*self.V[2]
                + self.g[0][g]*self.V[0] for g in glist}
            self.compSize = compH[glist[0]].shape[0]

            # Seed vector
            v0 = scipy.zeros(self.compSize)
            for i in range(min(10,len(v0))):
                v0[i] = 1.

            for g in glist:
                (self.eigenvalues[g][ren], eigenvectorstranspose) = \
                    scipy.sparse.linalg.eigsh(compH[g], neigs, v0=v0, which='SA', return_eigenvectors=True)
                    # scipy.sort(scipy.sparse.linalg.eigsh(compH[g], neigs, v0=v0,
                            # which='SA', return_eigenvectors=True))

                self.eigenvectors[g][ren] = eigenvectorstranspose.T

                if self.k == 1:
                    v = self.eigenvectors[g][ren][0]
                    v2 = (self.V[2]/self.L).dot(v)
                    norm = np.inner(v,v)
                    np.testing.assert_almost_equal(1,norm)
                    s = np.inner(v2,v)
                    self.vev[g][ren] = s

            return


        # Subset of low energy states
        subOpL = SubmatrixOperator.fromErange(self.basis,(0,ET))
        M = {g: self.h0 + g*self.V[4] + self.g[2][g]*self.V[2]
                + self.g[0][g]*self.V[0] for g in glist}
        Hraw = {g: subOpL.sub(M[g]) for g in glist}
        self.compSize = Hraw[glist[0]].shape[0]

        # Seed vector
        v0 = scipy.zeros(self.compSize)
        for i in range(min(10,len(v0))):
            v0[i] = 1.



        if ren=="rentails":
# We invert the matrices one by one to save memory
            DH2lL, DH2ll, DH3ll, DH2Ll =\
                self.computeDeltaH(ET=ET, EL=EL, ren=ren, eps=eps,
                    subbasisl=subbasisl, ELp=ELp, ELpp=ELpp, loc2=loc2,
                    loc3=loc3, loc3mix=loc3mix, nonloc3mix=nonloc3mix, memdbg=memdbg)

            # Saving memory
            if memsave:
                del self.basisH
                gc.collect()

            print("Diagonalizing matrices...")

            for g in glist:
                ML = DH2lL[g].todense()
                # For all tails
                MR = ML
                M = (DH2ll[g]-DH3ll[g]).todense()
                compH = Hraw[g].todense() + ML*scipy.linalg.solve(M, MR)

                self.eigenvalues[g][ren] = \
                    scipy.sort(scipy.sparse.linalg.eigsh(compH, neigs, v0=v0,
                        which='SA', return_eigenvectors=False))
                # self.eigenvectors[g][ren] = eigenvectorstranspose.T

                del ML, M, compH
                gc.collect()
            return

        else:
            raise ValueError()

        for g in glist:
            (self.eigenvalues[g][ren], eigenvectorstranspose) = \
                scipy.sparse.linalg.eigsh(compH[g], neigs, v0=v0,
                        which='SA', return_eigenvectors=True)

            self.eigenvectors[g][ren] = eigenvectorstranspose.T

            v = self.eigenvectors[g][ren][0]

            # if self.k == 1:
                # v = self.eigenvectors[g][ren][0]
                # v2 = (self.V[2]/self.L).dot(v)
                # s = np.inner(v2,v)
                # self.vev[g][ren] = s

