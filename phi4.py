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
from matrix import Matrix
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
        self.VHl = {}
        self.Vhl = {2:{}, 4:{}}
        self.Vll = {}
        self.V0V4 = {}
        self.V2V4 = {}
        self.V4V4 = {}
        # self.VhhDiagList = {}
        self.basisH = {}
        self.basish = {}
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


    def buildMatrix(self, Vlist, basis, lookupbasis, statePos=None,
            ignKeyErr=False, idxList=None, sumTranspose=False, subDiag=True):


        if basis.helper.nmax > lookupbasis.helper.nmax:
            helper = basis.helper
        else:
            helper = lookupbasis.helper

# Dictionary of positions of states
        # Contains also the P-reversed states
        # NOTE: using arrays is much less efficient!
        statePos = {}
        for i,state in enumerate(lookupbasis.stateList):
            statePos[tuple(helper.torepr2(state))] = i
            statePos[tuple(helper.torepr2(state)[::-1])] = i

        if idxList==None:
            idxList = range(basis.size)


        # Will construct the sparse matrix in the COO format and then convert it to CSC
        data = []
        row = []
        col = []

        for V in Vlist:
            for i in idxList:
                colpart, datapart = \
                    V.computeMatrixElements(basis,i,lookupbasis,
                            statePos=statePos, helper=helper, ignKeyErr=ignKeyErr)
                data += datapart
                col += colpart
                row += [i]*len(colpart)

        V = scipy.sparse.coo_matrix((data,(row,col)),
                shape=(basis.size,lookupbasis.size))

        if sumTranspose:
            # Add the matrix to its transpose and subtract the diagonal
            diag_V = scipy.sparse.spdiags(V.diagonal(),0,basis.size,basis.size).tocsc()
            if subDiag:
                return (V+V.transpose()-diag_V)
            else:
                return (V+V.transpose())
        else:
            return V



    def computePotential(self, k):
        """
        Builds the potential matrices and the free Hamiltonian. In the low-energy sector
        """

        basis = self.basis[k]

        self.h0[k] = Matrix(basis, basis,
                scipy.sparse.spdiags(basis.energyList,0,basis.size,basis.size))

        Vlist = {2:V2OpsHalf(basis), 4:V4OpsHalf(basis)}
        for n in (2,4):
            self.V[k][n] = Matrix(basis,basis,
                    self.buildMatrix(Vlist[n], basis, basis, sumTranspose=True)*self.L)

        # Construct the identity potential matrix
        idM = scipy.sparse.eye(basis.size)
        self.V[k][0] = Matrix(basis, basis, idM)*self.L


    def genHEBases(self, k, basisl, EL, ELpp):
        """ Generate a high-energy basis from a set of tails
        k: parity quantum number
        basisl: Basis instance containing the set of tails
        EL: maximal energy of the generated basis
        """

        self.basisl[k] = basisl

        for Emax,basis in zip((EL, ELpp),(self.basisH, self.basish)):
            # Generate all the operators between the selected states and the states
            # in the range [0, EL]

            if Emax == None:
                basis[k] = None
                continue

            Vlist = V4OpsSelectedFull(basisl, Emax)
            vectorset = set()

            for V in Vlist:
                for v in V.genBasis(basisl, Emax):
# Don't add twice states connected by parity inversion
                    if v not in vectorset and v[::-1] not in vectorset:
                        vectorset.add(v)

            helper = Vlist[0].helper

            # Basis of selected states with energy <= Emax
            basis[k] = Basis(k, [helper.torepr1(v) for v in vectorset], helper)


    # @profile
    def computeHEVs(self, k):
        """
        Compute the matrices involving the high-energy states below EL
        k: parity quantum number
        EL: local cutoff for DH2
        """

# NOTE matrix subscript notation:
# "l": selected low-energy state with energy <= ET
# "L": generic low-energy state with energy <= ET
# "h": selected high-energy state with energy <= EL'
# "H": selected high-energy states with energy <= EL

        #################################
        # Generate the VlH matrix
        #################################
        if self.basisH[k] != None:
            basis = self.basisl[k]
            lookupbasis = self.basisH[k]

            Vlist = V4OpsSelectedFull(basis, lookupbasis.Emax)

            self.VHl[k] = Matrix(basis, lookupbasis,
                    self.buildMatrix(Vlist, basis, lookupbasis)*self.L)


            ##############################
            # Generate the VHL matrix
            ##############################

            print("Computing VhL...")

            basis = self.basisH[k]
            lookupbasis = self.basis[k]

            Vlist = V4OpsSelectedFull(basis, lookupbasis.Emax)

            self.VLH[k] = Matrix(basis, lookupbasis,
                    self.buildMatrix(Vlist, basis, lookupbasis)*self.L)

            print("self.VLH[k] size", msize(self.VLH[k].M))

#################################
# Generate the local VhL matrices
#################################

        if self.basish[k] != None:
            basis = self.basisl[k]
            lookupbasis = self.basish[k]

            Vlist = V2OpsSelectedFull(basis, lookupbasis.Emax)
            self.Vhl[2][k] = Matrix(basis, lookupbasis,
                    self.buildMatrix(Vlist, basis, lookupbasis)*self.L)

            Vlist = V4OpsSelectedFull(basis, lookupbasis.Emax)
            self.Vhl[4][k] = Matrix(basis, lookupbasis,
                    self.buildMatrix(Vlist, basis, lookupbasis)*self.L)


    def computeLEVs(self, k):

        ###################################
        # Generate all the "local" matrices on the selected
        # low-energy states
        ###################################
        basis = self.basisl[k]

        self.Vll[k] = {}
        for n in (0,2,4):
            self.Vll[k][n] = self.V[k][n].sub(basis, basis)

        # Vlist = V6OpsSelectedHalf(basis)
        Vlist = V6OpsSelectedFull(basis, self.basis[k].Emax)

        self.Vll[k][6] = Matrix(basis,basis,
                self.buildMatrix(Vlist, basis, basis, ignKeyErr=True,
                    sumTranspose=False)*self.L)
        # print(self.Vll[k][6].M)

        ###################################
        # Generate all the "bilocal" matrices on the selected
        # low-energy states
        ###################################
        self.V0V4[k] = self.Vll[k][4]*self.L

        # Vlist = V2V4Ops1(basis)
        # self.V2V4[k] = Matrix(basis,basis,
                # self.buildMatrix(Vlist,basis,basis,ignKeyErr=True,
                    # sumTranspose=True, subDiag=False)*self.L**2)
        # Vlist = V2V4Ops2(basis)
        # self.V2V4[k] += Matrix(basis,basis,
                # self.buildMatrix(Vlist,basis,basis,ignKeyErr=True,
                        # sumTranspose=False)*self.L**2)

        print("Computing V2V4")

        Vlist = V2V4Ops(basis)
        self.V2V4[k] = Matrix(basis,basis,
                self.buildMatrix(Vlist,basis,basis,ignKeyErr=True,
                        sumTranspose=False)*self.L**2)



        # Vlist = V4V4Ops1(basis)
        # self.V4V4[k] = Matrix(basis,basis,
                # self.buildMatrix(Vlist,basis,basis,ignKeyErr=True,
                    # sumTranspose=True, subDiag=False)*self.L**2)
        # Vlist = V4V4Ops2(basis)
        # self.V4V4[k] += Matrix(basis,basis,
                # self.buildMatrix(Vlist,basis,basis,ignKeyErr=True,
                        # sumTranspose=False)*self.L**2)


        print("Computing V4V4")

        Vlist = V4V4Ops(basis)
        self.V4V4[k] = Matrix(basis,basis,
                self.buildMatrix(Vlist,basis,basis,ignKeyErr=True,
                        sumTranspose=False)*self.L**2)



    def computeDeltaH(self, k, ren, ET, eps, loc2=True, loc3=True, loc3mix=True,
            nonloc3mix=True, EL=None, ELp=None, ELpp=None, maxntails=None):
        """
        Compute the full DeltaH = DH2 * (DH2-DH3)^-1 * DH2  matrix
        """

        helper = self.basis[k].helper
        subbasisL = self.basis[k].sub(lambda v: helper.energy(v)<=ET)

        if ren=="rentails":
            basisl = self.basisl[k]
            helper = self.basisl[k].helper

            # Subset of the selected low energy states
            vectorlist = [v for v in basisl if helper.energy(v)<=ET]
            if maxntails != None:
                vectorlist = vectorlist[:maxntails]
            subbasisl = Basis(k, vectorlist, helper, orderEnergy=False)
            self.ntails = subbasisl.size

            DH2lL = self.computeDH2(subbasisl, k, ET=ET, EL=EL, eps=eps, loc2=loc2)
            DH2Ll = DH2lL.transpose()
            DH2ll = Matrix(subbasisl, subbasisL, DH2Ll).sub(subbasisl, subbasisl).M.tocsc()


            if self.DH3ll[k] == None:
                DH3ll = self.computeDH3(subbasisl, k, ET=ET, ELp=ELp, ELpp=ELpp, eps=eps,
                        loc3=loc3, loc3mix=loc3mix, nonloc3mix=nonloc3mix).M
            else:
                DH3ll = self.DH3ll[k].sub(subbasisl, subbasisl).M


            invM = scipy.sparse.linalg.inv(DH2ll-DH3ll)

            return DH2lL, invM, DH2Ll


        elif ren=="renloc":
            VV2 = renorm.renVV2(g4=self.g4, EL=ET, eps=eps).VV2

            # Dictionary of local renormalization coefficients for the g^2 term
            VLL = {}
            for n in (0,2,4):
                VLL[n] = self.V[k][n].sub(subbasisL, subbasisL)

            return (VV2[0]*VLL[0] + VV2[2]*VLL[2] + VV2[4]*VLL[4]).M



    def computeDH2(self, subbasisl, k, ET, EL, eps, loc2):

        helper = self.basis[k].helper
        # Subset of the full low energy states
        subbasisL = self.basis[k].sub(lambda v: helper.energy(v)<=ET)

        # Subsets of the selected high energy states
        helper = self.basisH[k].helper
        subbasisH = self.basisH[k].sub(lambda v: ET<helper.energy(v)<=EL)

        # Propagators on the high-energy states
        energyArr = array(subbasisH.energyList)
        propagatorH = scipy.sparse.spdiags(1/(eps-energyArr),0,subbasisH.size,
                        subbasisH.size)


        #############################
        # Construct the matrices needed to compute DeltaH
        #############################
        VHl = self.VHl[k].sub(subbasisl, subbasisH).M
        VlH = VHl.transpose()
        VLH = self.VLH[k].sub(subbasisH, subbasisL).M
        VHL = VLH.transpose()


        VlL = {}
        VLl = {}
        for n in (0,2,4):
            VlL[n] = self.V[k][n].sub(subbasisL, subbasisl).M
            VLl[n] = VlL[n].transpose()



        #######################################
        # Construct DH2
        #######################################
        # Dictionary of local renormalization coefficients for the g^2 term
        VV2 = renorm.renVV2(g4=self.g4, EL=EL, eps=eps).VV2

        # XXX Sorry, for now the subscripts are confusing, need to sort this out
        DH2lL = VHL*propagatorH*VlH*self.g4**2

        if loc2:
            DH2lL += VV2[0]*VlL[0] + VV2[2]*VlL[2] + VV2[4]*VlL[4]

        return DH2lL


    def precomputeDH3(self, k, ET, ELp, ELpp, eps, loc3=True, loc3mix=True, nonloc3mix=True):

        self.DH3ll[k] = self.computeDH3(self.basisl[k], k, ET, ELp, ELpp, eps, loc3, loc3mix, nonloc3mix)



    def computeDH3(self, subbasisl, k, ET, ELp, ELpp, eps, loc3, loc3mix, nonloc3mix):

        print("Computing DH3")

        helper = self.basish[k].helper
        subbasish = self.basish[k].sub(lambda v: ET<helper.energy(v)<=ELp)
        energyArr = array(subbasish.energyList)
        propagatorh = scipy.sparse.spdiags(1/(eps-energyArr),0,subbasish.size,
                subbasish.size)

# Local renormalization matrices used to approximate the sum over states with energy
# above EL
        Vhl = {}
        Vlh = {}
        for n in (2,4):
            Vhl[n] = self.Vhl[n][k].sub(subbasisl, subbasish).M
            Vlh[n] = Vhl[n].transpose()


# Low-energy matrix to fill in
        DH3ll = scipy.sparse.csc_matrix((subbasisl.size, subbasisl.size))


#########################################################################
# Add the "symmetric" contributions to DH3 by integrating out states with
# energy ET < E < ELp
#########################################################################
        basis = subbasish
        Vhhlist = V4OpsSelectedHalf(basis)


        chunklen = min(self.chunklen, basis.size)

        idxLists = [range(basis.size)[x:x+chunklen] for x in
                range(0, basis.size, chunklen)]


        for i,idxList in enumerate(idxLists):
            ##############################
            # Generate the Vhh matrix
            ##############################
            print("Doing chunk", i, "for Vhh")


# NOTE Trick to save memory: we never compute explicitly the full matrix Vhh
# TODO Do this more elegantly with iterators ?
            VhhHalfPart =  self.buildMatrix(Vhhlist, basis, basis,
                    ignKeyErr=True, sumTranspose=False, idxList=idxList)*self.L

            VhhDiagPart = scipy.sparse.spdiags(VhhHalfPart.diagonal(),0,basis.size,
                    basis.size)

            DH3llPart = Vhl[4]*propagatorh*VhhHalfPart*propagatorh*Vlh[4]*self.g4**3
            DH3llPart += DH3llPart.transpose()
            DH3llPart -= Vhl[4]*propagatorh*VhhDiagPart*propagatorh*Vlh[4]*self.g4**3
            DH3ll += DH3llPart

            del VhhHalfPart
            gc.collect()

#########################################################################
# Add the "mixed" contributions to DH3 by integrating out states with
# energy ET < E < ELp on one side and ELp < E < EL on the other
#########################################################################

        if nonloc3mix:
            # Subsets of the selected high energy states
            helper = self.basish[k].helper
            subbasisH = self.basish[k].sub(lambda v: ELp<helper.energy(v)<=ELpp)

            # Propagators on the high-energy states
            energyArr = array(subbasisH.energyList)
            propagatorH = scipy.sparse.spdiags(1/(eps-energyArr),0,subbasisH.size,
                    subbasisH.size)


            VHl = self.Vhl[4][k].sub(subbasisl, subbasisH).M
            VlH = VHl.transpose()


            basis = subbasish
            lookupbasis = subbasisH

            VhHlist = V4OpsSelectedFull(subbasish, ELpp)


            chunklen = min(self.chunklen, basis.size)

            idxLists = [range(basis.size)[x:x+chunklen] for x in
                range(0, basis.size, chunklen)]

            for i, idxList in enumerate(idxLists):

                print("doing chunk", i, "for VhH")

                VhHPart =  self.buildMatrix(VhHlist, basis, lookupbasis,
                        ignKeyErr=True, sumTranspose=False, idxList=idxList)*self.L

                DH3llPart = Vhl[4]*propagatorh*VhHPart*propagatorH*VlH*self.g4**3
                DH3llPart += DH3llPart.transpose()
                DH3ll += DH3llPart

                del VhHPart
                gc.collect()


#######################################################################
# Add the "mixed" local-nonlocal contributions to DH3 where some states
# are integrated from ET to ELp while others from EL to Infinity
########################################################################

        if loc3mix:
            VV2 = renorm.renVV2(g4=self.g4, EL=ELpp, eps=eps).VV2

            DH3llPart = Vhl[2]*VV2[2]*propagatorh*Vlh[4]*self.g4
            DH3llPart += Vhl[4]*VV2[4]*propagatorh*Vlh[4]*self.g4
            DH3llPart += DH3llPart.transpose()
            DH3ll += DH3llPart

#####################################################
# Add the "symmetric local" parts to DH3
#####################################################

        if loc3:
            Vll = {}
            for n in (0,2,4,6):
                Vll[n] = self.Vll[k][n].sub(subbasisl, subbasisl).M


            V0V4 = self.V0V4[k].sub(subbasisl,subbasisl).M
            V2V4 = self.V2V4[k].sub(subbasisl,subbasisl).M
            V4V4 = self.V4V4[k].sub(subbasisl,subbasisl).M


            DH3ll += V0V4*self.VV3.V0V4[ELp]*self.g4**3
            DH3ll += V2V4*self.VV3.V2V4[ELp]*self.g4**3
            DH3ll += V4V4*self.VV3.V4V4[ELp]*self.g4**3


            DH3ll += Vll[0]*self.VV3.VV3loc[0][ELp]*self.g4**3
            DH3ll += Vll[2]*self.VV3.VV3loc[2][ELp]*self.g4**3
            DH3ll += Vll[4]*self.VV3.VV3loc[4][ELp]*self.g4**3
            DH3ll += Vll[6]*self.VV3.VV3loc[6][ELp]*self.g4**3

        return Matrix(subbasisl, subbasisl, DH3ll)



    def calcVV3(self, ELplist, eps):
        print("Calculating VVV renorm coefficients")
        self.VV3 =  renorm.renVV3(m=self.m, eps=eps, ETlist=ELplist)


        # print("Bilocal g^3 ren coefficients: ", self.VV3.V0V4, self.VV3.V2V4, self.VV3.V4V4)

        # print("Local g^3 ren coefficients: ", self.VV3.VV3loc)


    def setCouplings(self, g0=0, g2=0, g4=0):
        self.g0 = g0
        self.g2 = g2
        self.g4 = g4



    def computeEigval(self, k, ET, ren, EL=None, ELp=None, ELpp=None, loc2=True,
            eps=None, neigs=10, maxntails=None, loc3=True, loc3mix=True, nonloc3mix=True):
        """ Compute the eigenvalues for sharp cutoff ET and local cutoff EL
        k: parity quantum number
        ET: ET
        EL: EL. Should be set to None for "raw" and to EL for "renlocal"
        ren: type of eigenvalue. Can be "raw", "renlocal" or "rentails"
        eps: epsilon parameter in the propagator
        neigs: number of eigenvalues to compute
        """

        # Subset of low energy states
        helper = self.basis[k].helper

        subbasisL = self.basis[k].sub(lambda v: helper.energy(v)<=ET)
        self.subbasisL = subbasisL

        Hraw = (self.h0[k] + self.g4*self.V[k][4]).sub(subbasisL, subbasisL).M

        if ren=="raw":
            compH = Hraw

        elif ren=="rentails":
            DH2lL, invM, DH2Ll = self.computeDeltaH(k=k, ET=ET, EL=EL, ren=ren, eps=eps,
                    maxntails=maxntails, ELp=ELp, ELpp=ELpp, loc2=loc2,
                    loc3=loc3, loc3mix=loc3mix, nonloc3mix=nonloc3mix)
            compH = LinearOperator(Hraw.shape, lambda v: Hraw*v + DH2lL*invM*DH2Ll*v)

        elif ren=="renloc":
            DeltaH = self.computeDeltaH(k=k, ET=ET, EL=EL, ren=ren, eps=eps,
                    maxntails=maxntails, ELp=ELp, ELpp=ELpp, loc2=loc2,
                    loc3=loc3, loc3mix=loc3mix, nonloc3mix=nonloc3mix)
            compH = (Hraw + DeltaH)
        else:
            raise ValueError()

        self.compSize = Hraw.shape[0]

        # Seed vector
        v0 = scipy.zeros(subbasisL.size)
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

