import scipy
import scipy.sparse.linalg
import scipy.sparse
import math
from math import factorial
from statefuncs import Basis, Helper
from oscillators import *
from collections import Counter
from operator import attrgetter
import renorm
import gc
from matrix import *
from scipy import exp, pi, array
from scipy.sparse.linalg import LinearOperator



def msize(m):
    form =  m.getformat()
    if form=="coo":
        return (m.col.nbytes+m.row.nbytes+m.data.nbytes)/1000000
    elif form == "csr" or form=="csc":
        return (m.data.nbytes+m.indptr.nbytes+m.indices.nbytes)/1000000
    else:
        raise ValueError(form+" not implemented")



class Phi4():
    """ main class """
    def __init__(self, m, L):

        self.m = m
        self.L = L
# Maximum dimension of the chunks for computing Vhh
        self.chunklen = 20000

        self.basis = {}
        self.h0 = {}
        self.V = {k:{} for k in {-1,1}}
        self.DeltaH = {}
        self.VLH = {}
        self.VHl = {k:{} for k in (-1,1)}
        self.Vll = {k:{} for k in (-1,1)}
        self.VlL = {k:{} for k in (-1,1)}
        self.V0V4 = {}
        self.V2V4 = {}
        self.V4V4 = {}
        # self.VhhDiagList = {}
        self.basisH = {}
        self.basisl = {}
        self.DH3ll = {k:None for k in (-1,1)}


# "raw" are the raw eigenvalues
# "renloc" are the eigenvalues where the DH2, DH3 corrections are only computed in
# the local approximation
# "rentails" are the eigenvalues where the DH2, DH3 corrections are computed exactly up
# to EL
        self.eigenvalues = {"raw":{}, "renloc":{}, "rentails":{}}
        self.eigenvectors = {"raw":{}, "renloc":{}, "rentails":{}}

        scipy.set_printoptions(precision=15)


    def buildBasis(self, Emax, occmax=None):
        """ Builds the full Hilbert space basis up to cutoff Emax """
        self.basis = Basis.fromScratch(m=self.m, L=self.L, Emax=Emax, occmax=occmax)


    def computePotential(self, k):
        """
        Builds the potential matrices and the free Hamiltonian. In the low-energy sector
        """

        basis = self.basis[k]

        self.h0[k] = scipy.sparse.spdiags(basis.energyList,0,basis.size,basis.size)

        c = MatrixConstructor(basis, basis)
        Vlist = {2:V2OpsHalf(basis), 4:V4OpsHalf(basis)}
        for n in (2,4):
            print("computing V", n)
            self.V[k][n] = c.buildMatrix(Vlist[n], sumTranspose=True)*self.L
        del c

        # Construct the identity potential matrix
        self.V[k][0] = scipy.sparse.eye(basis.size)*self.L


    # @profile
    def genHEBasis(self, k, basisl, EL, ELp, ELpp):
        """ Generate a high-energy basis from a set of tails
        k: parity quantum number
        basisl: Basis instance containing the set of tails
        EL: maximal energy of the generated basis for DH2
        ELpp: maximal energy of the generated basis for DH3
        """

        self.basisl[k] = basisl

        self.EL = EL
        self.ELp = ELp
        self.ELpp = ELpp
        Emax = max(EL, ELp, ELpp)

        # Generate all the operators between the selected states and the states
        # in the range [0, Emax]
        Vlist = V4OpsSelectedFull(basisl, Emax)
        vectorset = set()

        for V in Vlist:
            for v in V.genBasis(basisl, Emax):
                # Don't add twice states connected by parity inversion
                if v not in vectorset and v[::-1] not in vectorset:
                    vectorset.add(v)

        helper = Vlist[0].helper

        # Basis of selected states with energy <= Emax
        self.basisH[k] = Basis(k, (helper.torepr1(v) for v in vectorset), helper)


# XXX We could compute either Vhl or VHl to save time
    def computeHEVs(self, k):
        """
        Compute the matrices involving the high-energy states below EL
        k: parity quantum number
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

        basis = self.basisl[k]
        lookupbasis = self.basisH[k]
        Emax = lookupbasis.Emax

        c = MatrixConstructor(basis, lookupbasis)

        Vlist = V4OpsSelectedFull(basis, Emax)
        self.VHl[k][4] = c.buildMatrix(Vlist)*self.L

        Vlist = V2OpsSelectedFull(basis, Emax)
        self.VHl[k][2] = c.buildMatrix(Vlist)*self.L

        del c

        ##############################
        # Generate the VHL matrix
        ##############################

        print("Computing VLH...")

        basis = self.basisH[k]
        lookupbasis = self.basis[k]
        idxList = basis.irange((0, Emax))

        c = MatrixConstructor(basis, lookupbasis)

        Vlist = V4OpsSelectedFull(basis, lookupbasis.Emax, idxList=idxList)
        self.VLH[k] = c.buildMatrix(Vlist, idxList=idxList)*self.L

        del c

        print("self.VLH[k] size", msize(self.VLH[k]))


    def computeLEVs(self, k):

        ###################################
        # Generate all the "local" matrices on the selected
        # low-energy states
        ###################################

        subOp = SubmatrixOperator.fromSubbasis(self.basis[k], self.basisl[k])

        self.Vll[k] = {}
        for n in (0,2,4):
            self.VlL[k][n] = subOp.subcolumns(self.V[k][n])
            self.Vll[k][n] = subOp.subrows(self.VlL[k][n])


        basis = self.basisl[k]
        c = MatrixConstructor(basis, basis)

        Vlist = V6OpsSelectedFull(basis, basis.Emax)
        self.Vll[k][6] = c.buildMatrix(Vlist, ignKeyErr=True, sumTranspose=False)*self.L

        ###################################
        # Generate all the "bilocal" matrices on the selected
        # low-energy states
        ###################################
        self.V0V4[k] = self.Vll[k][4]*self.L

        print("Computing V2V4")

        Vlist = V2V4Ops(basis)
        self.V2V4[k] = c.buildMatrix(Vlist,ignKeyErr=True,sumTranspose=False)*self.L**2

        print("Computing V4V4")

        Vlist = V4V4Ops(basis)
        self.V4V4[k] = c.buildMatrix(Vlist,ignKeyErr=True,sumTranspose=False)*self.L**2

        del c


    def computeDeltaH(self, k, ren, ET, eps, loc2=True, loc3=True, loc3mix=True,
            nonloc3mix=True, EL=None, ELp=None, ELpp=None, subbasisl=None):
        """
        Compute the full DeltaH = DH2 * (DH2-DH3)^-1 * DH2  matrix
        subbasisl: explicit set of tails used in the computation. If not specified,
        the tails in self.basisl below ET are used
        """

        # These operators allow to take submatrices according to a given energy range
        subOPL = SubmatrixOperator.fromErange(self.basis[k], (0, ET))

        if ren=="rentails":
            if subbasisl==None:
                subOPl = SubmatrixOperator.fromErange(self.basisl[k], (0,ET))
            else:
                subOPl = SubmatrixOperator.fromSubbasis(self.basisl[k], subbasisl)

            self.ntails = len(subOPl.idxList)

            DH2ll, DH2Ll = self.computeDH2(k, ET=ET, EL=EL, eps=eps, loc2=loc2)

            # Extract the submatrix in the right energy range
            DH2Ll = subOPL.subcolumns(subOPl.subrows(DH2Ll))
            DH2lL = DH2Ll.transpose()

            DH2ll = subOPl.sub(DH2ll)

            DH3ll = self.computeDH3(k, ET=ET, ELp=ELp, ELpp=ELpp, eps=eps,
                        loc3=loc3, loc3mix=loc3mix, nonloc3mix=nonloc3mix)
            DH3ll = subOPl.sub(DH3ll)

            invM = scipy.sparse.linalg.inv((DH2ll-DH3ll).tocsc())

            return DH2lL, invM, DH2Ll


        elif ren=="renloc":

            # Dictionary of local renormalization coefficients for the g^2 term
            VV2 = renorm.renVV2(g4=self.g4, EL=ET, eps=eps).VV2

            VLL = {}
            for n in (0,2,4):
                VLL[n] = subOPL.sub(self.V[k][n])

            return VV2[0]*VLL[0] + VV2[2]*VLL[2] + VV2[4]*VLL[4]


    def computeDH2(self, k, ET, EL, eps, loc2):

        propagatorH = self.basisH[k].propagator(eps, ET, EL)

        VHl = self.VHl[k][4]
        VlH = VHl.transpose()
        VLH = self.VLH[k]

        VlL = {}
        VLl = {}
        Vll = {}
        for n in (0,2,4):
            VlL[n] = self.VlL[k][n]
            VLl[n] = VlL[n].transpose()
            Vll[n] = self.Vll[k][n]

        #######################################
        # Construct DH2
        #######################################
        # Dictionary of local renormalization coefficients for the g^2 term
        VV2 = renorm.renVV2(g4=self.g4, EL=EL, eps=eps).VV2

        DH2Ll = VHl*propagatorH*VLH*self.g4**2
        DH2ll = VHl*propagatorH*VlH*self.g4**2

        if loc2:
            DH2Ll += VV2[0]*VLl[0] + VV2[2]*VLl[2] + VV2[4]*VLl[4]
            DH2ll += VV2[0]*Vll[0] + VV2[2]*Vll[2] + VV2[4]*Vll[4]

        return DH2ll, DH2Ll


    def computeDH3(self, k, ET, ELp, ELpp, eps, loc3, loc3mix, nonloc3mix):

        print("Computing DH3")

        sizel = self.basisl[k].size
        DH3ll = scipy.sparse.csc_matrix((sizel, sizel))

        VHl = {}
        VlH = {}
        for n in (2,4):
            VHl[n] = self.VHl[k][n]
            VlH[n] = VHl[n].transpose()

        basis = self.basisH[k]
        # List of basis elements on which we will cycle
        fullIdxList = basis.irange((ET,ELp))
        idxlen = len(fullIdxList)
        chunklen = min(self.chunklen, idxlen)
        idxLists = [fullIdxList[x:x+chunklen] for x in range(0, idxlen, chunklen)]

        print("Total number of chunks: ", len(idxLists))

#########################################################################
# Add the "symmetric" contributions to DH3 by integrating out states with
# energy ET < E < ELp
#########################################################################

# Propagator and projector over the states between ET and ELp
        propagatorH = basis.propagator(eps, ET, ELp)

        Vlist = V4OpsSelectedHalf(basis, Emax=ELp, idxList=fullIdxList)

        c = MatrixConstructor(basis, basis, (ET, ELp))

        ##############################
        # Generate the Vhh matrix
        ##############################
        for i,idxList in enumerate(idxLists):
            print("Doing chunk", i, "for Vhh")

            VhhHalfPart =  c.buildMatrix(Vlist, ignKeyErr=True, sumTranspose=False,
                    idxList=idxList)*self.L

            VhhDiagPart = scipy.sparse.spdiags(VhhHalfPart.diagonal(),0,basis.size,
                    basis.size)

            DH3llPart = VHl[4]*propagatorH*VhhHalfPart*propagatorH*VlH[4]*self.g4**3
            DH3llPart += DH3llPart.transpose()
            DH3llPart -= VHl[4]*propagatorH*VhhDiagPart*propagatorH*VlH[4]*self.g4**3
            DH3ll += DH3llPart

            del VhhHalfPart

        del c

#########################################################################
# Add the "mixed" contributions to DH3 by integrating out states with
# energy ET < E < ELp on one side and ELp < E < EL on the other
#########################################################################

        if nonloc3mix:
            c = MatrixConstructor(basis, basis, (ELp, ELpp))

# Propagator and projector over the states between ELp and ELpp
            propagatorHH = basis.propagator(eps, ELp, ELpp)

            VhHlist = V4OpsSelectedFull(basis, ELpp, idxList=fullIdxList)

            for i, idxList in enumerate(idxLists):
                print("doing chunk", i, "for VhH")

                VhHPart = c.buildMatrix(VhHlist, ignKeyErr=True, sumTranspose=False,
                        idxList=idxList)*self.L

                DH3llPart = VHl[4]*propagatorH*VhHPart*propagatorHH*VlH[4]*self.g4**3
                DH3llPart += DH3llPart.transpose()
                DH3ll += DH3llPart

                del VhHPart
                gc.collect()

            del c


#######################################################################
# Add the "mixed" local-nonlocal contributions to DH3 where some states
# are integrated from ET to ELp while others from ELpp to Infinity
########################################################################

        if loc3mix:
            VV2 = renorm.renVV2(g4=self.g4, EL=ELpp, eps=eps).VV2

            DH3llPart = VHl[2]*VV2[2]*propagatorH*VlH[4]*self.g4
            DH3llPart += VHl[4]*VV2[4]*propagatorH*VlH[4]*self.g4
            DH3llPart += DH3llPart.transpose()
            DH3ll += DH3llPart

#####################################################
# Add the "symmetric local" parts to DH3
#####################################################

        if loc3:
            Vll = {}
            for n in (0,2,4,6):
                Vll[n] = self.Vll[k][n]

            V0V4 = self.V0V4[k]
            V2V4 = self.V2V4[k]
            V4V4 = self.V4V4[k]


            DH3ll += V0V4*self.VV3.V0V4[ELp]*self.g4**3
            DH3ll += V2V4*self.VV3.V2V4[ELp]*self.g4**3
            DH3ll += V4V4*self.VV3.V4V4[ELp]*self.g4**3


            DH3ll += Vll[0]*self.VV3.VV3loc[0][ELp]*self.g4**3
            DH3ll += Vll[2]*self.VV3.VV3loc[2][ELp]*self.g4**3
            DH3ll += Vll[4]*self.VV3.VV3loc[4][ELp]*self.g4**3
            DH3ll += Vll[6]*self.VV3.VV3loc[6][ELp]*self.g4**3

        return DH3ll


    def calcVV3(self, ELplist, eps):
        print("Calculating VVV renorm coefficients")
        self.VV3 =  renorm.renVV3(m=self.m, eps=eps, ETlist=ELplist)


    def setCouplings(self, g0=0, g2=0, g4=0):
        self.g0 = g0
        self.g2 = g2
        self.g4 = g4

    def computeEigval(self, k, ET, ren, EL=None, ELp=None, ELpp=None, loc2=True,
            eps=None, neigs=10, subbasisl=None, loc3=True, loc3mix=True, nonloc3mix=True):
        """ Compute the eigenvalues for sharp cutoff ET and local cutoff EL
        k: parity quantum number
        ET: ET
        EL: EL. Should be set to None for "raw" and to EL for "renlocal"
        ren: type of eigenvalue. Can be "raw", "renlocal" or "rentails"
        eps: epsilon parameter in the propagator
        neigs: number of eigenvalues to compute
        """

        # Subset of low energy states
        subOpL = SubmatrixOperator.fromErange(self.basis[k],(0,ET))
        M = self.h0[k] + self.g4*self.V[k][4]
        Hraw = subOpL.sub(M)
        self.compSize = Hraw.shape[0]

        if ren=="raw":
            compH = Hraw

        elif ren=="rentails":
            DH2lL, invM, DH2Ll = self.computeDeltaH(k=k, ET=ET, EL=EL, ren=ren, eps=eps,
                    subbasisl=subbasisl, ELp=ELp, ELpp=ELpp, loc2=loc2,
                    loc3=loc3, loc3mix=loc3mix, nonloc3mix=nonloc3mix)
            compH = LinearOperator(Hraw.shape, lambda v: Hraw*v + DH2lL*invM*DH2Ll*v)

        elif ren=="renloc":
            DeltaH = self.computeDeltaH(k=k, ET=ET, ren=ren, eps=eps, loc2=loc2)
            compH = (Hraw + DeltaH)

        else:
            raise ValueError()

        # Seed vector
        v0 = scipy.zeros(self.compSize)
        for i in range(min(10,len(v0))):
            v0[i] = 1.

        (self.eigenvalues[ren][k], eigenvectorstranspose) = \
            scipy.sparse.linalg.eigsh(compH, neigs, v0=v0,
                        which='SA', return_eigenvectors=True)

        self.eigenvectors[ren][k] = eigenvectorstranspose.T

    def vacuumE(self, ren):
        return self.eigenvalues[ren][1][0]
        # The vacuum is K-even

    def spectrum(self, k, ren):
        eigs = self.eigenvalues[ren][k]
        # Subtract vacuum energies
        if k==1:
            return scipy.array([x-self.vacuumE(ren=ren) for x in eigs[1:]])
        elif k==-1:
            return scipy.array([x-self.vacuumE(ren=ren) for x in eigs])

