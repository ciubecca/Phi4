import scipy
from profile_support import *
import scipy.sparse.linalg
import scipy.sparse
import math
from math import factorial
from statefuncs import Basis, Helper
from oscillators import *
from bilocal import *
from collections import Counter
from operator import attrgetter
import renorm
import gc
from matrix import *
from scipy import exp, pi, array
from scipy.sparse.linalg import LinearOperator
from sys import getsizeof as sizeof
from multiprocessing import Pool
from functools import partial
import psutil


def msize(m):
    form =  m.getformat()
    if form=="coo":
        return (m.col.nbytes+m.row.nbytes+m.data.nbytes)/1000000
    elif form == "csr" or form=="csc":
        return (m.data.nbytes+m.indptr.nbytes+m.indices.nbytes)/1000000
    else:
        raise ValueError(form+" not implemented")

def printMemoryUsage():
    print("Memory of parent process in MB: ", memory_usage()[0])
    print("Parent + children  memory in MB: ",
            memory_usage()[0] +
            sum(memory_usage(child.pid)[0] for child in psutil.Process().children()))

class Phi4():
    """ main class """
    def __init__(self, m, L):

        self.m = m
        self.L = L
# Maximum dimension of the chunks for computing Vhh
        self.chunklen = 20000
# Number of processes to spawn
        self.nproc = 4
        print("Number of parallel processes:", self.nproc)
# Initialize the child processes. NOTE Doing it now could save memory
        self.pool = Pool(self.nproc)

        self.basis = {}
        self.h0 = {}
        self.V = {k:{} for k in {-1,1}}
        self.DeltaH = {}
        self.VLH = {}
        self.VHL = {}
        self.VHl = {k:{} for k in (-1,1)}
        self.Vll = {k:{} for k in (-1,1)}
        self.VlL = {k:{} for k in (-1,1)}
        self.V0V4 = {}
        self.V2V4 = {}
        self.V4V4 = {}
        self.basisH = {}
        self.basisl = {}
        self.DH3ll = {k:None for k in (-1,1)}


# "raw" are the raw eigenvalues
# "renloc" are the eigenvalues where the DH2, DH3 corrections are only computed in
# the local approximation
# "rentails" are the eigenvalues where the DH2, DH3 corrections are computed exactly up
# to EL

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
            self.V[k][n] = c.buildMatrix(Vlist[n], sumTranspose=True)*self.L
        del c

        # Construct the identity potential matrix
        self.V[k][0] = scipy.sparse.eye(basis.size)*self.L


    # @profile
    def genHEBasis(self, k, EL, ELp, ELpp):
        """ Generate a high-energy basis from a set of tails
        k: parity quantum number
        basisl: Basis instance containing the set of tails
        EL: maximal energy of the generated basis for DH2
        ELpp: maximal energy of the generated basis for DH3
        """

        basisl = self.basisl[k]

        self.EL = EL
        self.ELp = ELp
        self.ELpp = ELpp
        Emax = max(EL, ELp, ELpp)

        # Generate all the operators between the selected states and the states
        # in the range [0, Emax]
        Vlist = V4OpsSelectedFull(basisl, Emax)
        vectorset = set()

        for V in Vlist:
            for v in V.yieldBasis(basisl, Emax):
                # Don't add twice states connected by parity inversion
                if v not in vectorset and v[::-1] not in vectorset:
                    vectorset.add(v)

        helper = Vlist[0].helper

        # Basis of selected states with energy <= Emax
        self.basisH[k] = Basis(k, vectorset, helper, repr1=False, repr1Emax=ELp)
        # self.basisH[k] = Basis(k, map(helper.torepr1, vectorset), helper)


# XXX We could compute either Vhl or VHl to save time
    # @profile
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

        print("Computing VHL...")

        basis = self.basis[k]
        lookupbasis = self.basisH[k]

        # We only need this matrix for DH2, not for DH3
        Emax = self.EL
        Erange = (0,Emax)

        c = MatrixConstructor(basis, lookupbasis, Erange=Erange)

        # We just need the full operator
        Vlist = V4OpsSelectedFull(basis, Emax)
        self.VHL[k] = c.buildMatrix(Vlist, ignKeyErr=True)*self.L
        self.VLH[k] = self.VHL[k].transpose()

        del c



    # @profile
    def computeLEVs(self, k, basisl, loc3=True):

        ###################################
        # Generate all the "local" matrices on the selected
        # low-energy states
        ###################################

        self.basisl[k] = basisl

        subOp = SubmatrixOperator.fromSubbasis(self.basis[k], self.basisl[k])

        self.Vll[k] = {}
        for n in (0,2,4):
            self.VlL[k][n] = subOp.subcolumns(self.V[k][n])
            self.Vll[k][n] = subOp.subrows(self.VlL[k][n])


        basis = self.basisl[k]
        c = MatrixConstructor(basis, basis)

        Vlist = V6OpsSelectedFull(basis, basis.Emax)
        self.Vll[k][6] = c.buildMatrix(Vlist, ignKeyErr=True,
                sumTranspose=False)*self.L

        ###################################
        # Generate all the "bilocal" matrices on the selected
        # low-energy states
        ###################################

        if loc3:

            basis.helper.calcOscEnergyDict()

            self.V0V4[k] = self.Vll[k][4]*self.L

            print("Computing V2V4")

            Vlist = V2V4Ops(basis)
            self.V2V4[k] = c.buildMatrix(Vlist,ignKeyErr=True,
                    sumTranspose=False)*self.L**2

            print("Computing V4V4")

            sizel = self.basisl[k].size
            self.V4V4[k] = scipy.sparse.csc_matrix((sizel, sizel))

            for nd1 in (0,1,2,3,4):
                for nd2 in (0,1,2,3,4):
                    Vlist = V4V4Ops(basis,nd1,nd2)
                    self.V4V4[k] += c.buildMatrix(Vlist,ignKeyErr=True,
                        sumTranspose=False)*self.L**2

            del c

    # @profile
    def computeDeltaH(self, k, ren, ET, eps, loc2=True, loc3=True, loc3mix=True,
            nonloc3mix=True, EL=None, ELp=None, ELpp=None, subbasisl=None,
            memdbg=False):
        """
        Compute the full DeltaH = DH2 * (DH2-DH3)^-1 * DH2  matrix
        subbasisl: explicit set of tails used in the computation. If not specified,
        the tails in self.basisl below ET are used
        """

        glist = self.glist

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
            DH2Ll = {g: subOPL.subcolumns(subOPl.subrows(DH2Ll[g])) for g in glist}
            DH2lL = {g: DH2Ll[g].transpose() for g in glist}

            DH2ll = {g: subOPl.sub(DH2ll[g]) for g in glist}

            DH3ll = self.computeDH3(k, ET=ET, ELp=ELp, ELpp=ELpp, eps=eps,
                        loc3=loc3, loc3mix=loc3mix, nonloc3mix=nonloc3mix,
                        memdbg=memdbg)
            DH3ll = {g: subOPl.sub(DH3ll[g]) for g in glist}

            # XXX This should not be needed anymore
            del self.basisH[k]
            gc.collect()

            # XXX Note: this could take a lot of memory if there are many tails,
            # because the inverse matrix is not sparse
            if memdbg: printMemoryUsage()
            # Convert to dense, because the inverse will be dense
            #TODO Make this work
            # invM = scipy.linalg.inv((DH2ll-DH3ll).todense())

            return DH2lL, DH2ll, DH3ll, DH2Ll


        elif ren=="renloc":

            ret = {}

            VLL = {}
            for n in (0,2,4):
                VLL[n] = subOPL.sub(self.V[k][n])

            for g in glist:
                # Dictionary of local renormalization coefficients for the g^2 term
                VV2 = renorm.renVV2(g4=g, EL=ET, eps=eps[g]).VV2

                ret[g] = VV2[0]*VLL[0] + VV2[2]*VLL[2] + VV2[4]*VLL[4]

            return ret


    def computeDH2(self, k, ET, EL, eps, loc2):

        glist=self.glist

        propagatorH = {g: self.basisH[k].propagator(eps[g], ET, EL) for g in glist}

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

        DH2Ll = {}
        DH2ll = {}

        for g in glist:
            VV2 = renorm.renVV2(g4=g, EL=EL, eps=eps[g]).VV2

            DH2Ll[g] = VHl*propagatorH[g]*VLH*g**2
            DH2ll[g] = VHl*propagatorH[g]*VlH*g**2

            if loc2:
                DH2Ll[g] += VV2[0]*VLl[0] + VV2[2]*VLl[2] + VV2[4]*VLl[4]
                DH2ll[g] += VV2[0]*Vll[0] + VV2[2]*Vll[2] + VV2[4]*Vll[4]

        return DH2ll, DH2Ll

    # @profile
    def computeDH3(self, k, ET, ELp, ELpp, eps, loc3, loc3mix, nonloc3mix, memdbg):

        glist = self.glist

        print("Computing DH3")

        sizel = self.basisl[k].size
        DH3ll = {g: scipy.sparse.csc_matrix((sizel, sizel)) for g in glist}

        if memdbg: printMemoryUsage()

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
        propagatorh = {g: basis.propagator(eps[g], ET, ELp) for g in glist}

        # XXX Move this inside the cycle and restrict idxList?
        Vhhlist = V4OpsSelectedHalf(basis, Emax=ELp, idxList=fullIdxList)

        c = MatrixConstructor(basis, basis, (ET, ELp))

        ##############################
        # Generate the Vhh matrix
        ##############################

        if self.nproc>1:
            f = partial(buildMatrixChunk, L=self.L, glist=glist, Vlist=Vhhlist,
                    VL=VHl[4], VR=VlH[4], propL=propagatorh, propR=propagatorh,
                    c=c, subDiag=True)

            DH3llchunks = self.pool.map(f, idxLists)

            for g in glist:
                DH3ll[g] += sum(chunk[g] for chunk in DH3llchunks)

        else:
            for i,idxList in enumerate(idxLists):
                print("Doing chunk", i, "for Vhh")

                VhhHalfPart =  c.buildMatrix(Vhhlist, ignKeyErr=True,
                        sumTranspose=False, idxList=idxList)*self.L

                VhhDiagPart = scipy.sparse.spdiags(VhhHalfPart.diagonal(),0,basis.size,
                        basis.size)

                for g in glist:
                    DH3llPart = VHl[4]*propagatorh[g]*VhhHalfPart*propagatorh[g]\
                        *VlH[4]*g**3
                    DH3llPart += DH3llPart.transpose()
                    DH3llPart -= VHl[4]*propagatorh[g]*VhhDiagPart*propagatorh[g]\
                        *VlH[4]*g**3
                    DH3ll[g] += DH3llPart

                del VhhHalfPart

        del Vhhlist
        del c

#########################################################################
# Add the "mixed" contributions to DH3 by integrating out states with
# energy ET < E < ELp on one side and ELp < E < EL on the other
#########################################################################
        if memdbg: printMemoryUsage()

        if nonloc3mix:
            c = MatrixConstructor(basis, basis, (ELp, ELpp))

# Propagator and projector over the states between ELp and ELpp
            propagatorH = {g: basis.propagator(eps[g], ELp, ELpp) for g in glist}

            # XXX Move this inside the cycle and restrict idxList?
            VHhlist = V4OpsSelectedFull(basis, ELpp, idxList=fullIdxList)

            if self.nproc>1:
                f = partial(buildMatrixChunk, L=self.L, glist=glist, Vlist=VHhlist,
                        VL=VHl[4], VR=VlH[4], propL=propagatorh, propR=propagatorH,
                        c=c, subDiag=False)

                DH3llchunks = self.pool.map(f, idxLists)

                for g in glist:
                    DH3ll[g] += sum(chunk[g] for chunk in DH3llchunks)

            else:
                for i, idxList in enumerate(idxLists):
                    print("doing chunk", i, "for VhH")

                    VHhPart = c.buildMatrix(VHhlist, ignKeyErr=True, sumTranspose=False,
                            idxList=idxList)*self.L

                    for g in glist:
                        DH3llPart = VHl[4]*propagatorh[g]*VHhPart*propagatorH[g]*VlH[4]\
                            *g**3
                        DH3llPart += DH3llPart.transpose()
                        DH3ll[g] += DH3llPart

                    del VHhPart

            del c
            del VHhlist

#######################################################################
# Add the "mixed" local-nonlocal contributions to DH3 where some states
# are integrated from ET to ELp while others from ELpp to Infinity
########################################################################

        if loc3mix:
            for g in glist:
                VV2 = renorm.renVV2(g4=g, EL=ELpp, eps=eps[g]).VV2

                DH3llPart = VHl[2]*VV2[2]*propagatorh[g]*VlH[4]*g
                DH3llPart += VHl[4]*VV2[4]*propagatorh[g]*VlH[4]*g
                DH3llPart += DH3llPart.transpose()
                DH3ll[g] += DH3llPart

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
        self.glist = glist

        self.eigenvalues = {g: {"raw":{}, "renloc":{}, "rentails":{}} for g in glist}
        self.eigenvectors = {g: {"raw":{}, "renloc":{}, "rentails":{}} for g in glist}


    def computeEigval(self, k, ET, ren, EL=None, ELp=None, ELpp=None, loc2=True,
            eps=None, neigs=10, subbasisl=None, loc3=True, loc3mix=True,
            nonloc3mix=True, memdbg=False):
        """ Compute the eigenvalues for sharp cutoff ET and local cutoff EL
        k: parity quantum number
        ET: ET
        EL: EL. Should be set to None for "raw" and to EL for "renlocal"
        ren: type of eigenvalue. Can be "raw", "renlocal" or "rentails"
        eps: epsilon parameter in the propagator
        neigs: number of eigenvalues to compute
        """

        glist = self.glist

        # Subset of low energy states
        subOpL = SubmatrixOperator.fromErange(self.basis[k],(0,ET))
        M = {g: self.h0[k] + g*self.V[k][4] for g in glist}
        Hraw = {g: subOpL.sub(M[g]) for g in glist}
        self.compSize = Hraw[glist[0]].shape[0]

        # Seed vector
        v0 = scipy.zeros(self.compSize)
        for i in range(min(10,len(v0))):
            v0[i] = 1.

        if ren=="raw":
            compH = Hraw

        elif ren=="renloc":
            DeltaH = self.computeDeltaH(k=k, ET=ET, ren=ren, eps=eps, loc2=loc2)
            compH = {g: (Hraw[g] + DeltaH[g]) for g in glist}

        elif ren=="rentails":
# We invert the matrices one by one to save memory
            DH2lL, DH2ll, DH3ll, DH2Ll =\
                self.computeDeltaH(k=k, ET=ET, EL=EL, ren=ren, eps=eps,
                    subbasisl=subbasisl, ELp=ELp, ELpp=ELpp, loc2=loc2,
                    loc3=loc3, loc3mix=loc3mix, nonloc3mix=nonloc3mix, memdbg=memdbg)

            for g in glist:
                invM = scipy.sparse.linalg.inv((DH2ll[g]-DH3ll[g]).tocsc())

                compH = LinearOperator(Hraw[g].shape,
                        lambda v: Hraw[g]*v + DH2lL[g]*invM*DH2Ll[g]*v)

                (self.eigenvalues[g][ren][k], eigenvectorstranspose) = \
                    scipy.sparse.linalg.eigsh(compH, neigs, v0=v0,
                        which='SA', return_eigenvectors=True)

                self.eigenvectors[g][ren][k] = eigenvectorstranspose.T

            return

        else:
            raise ValueError()

        for g in glist:

            (self.eigenvalues[g][ren][k], eigenvectorstranspose) = \
                scipy.sparse.linalg.eigsh(compH[g], neigs, v0=v0,
                        which='SA', return_eigenvectors=True)

            self.eigenvectors[g][ren][k] = eigenvectorstranspose.T

