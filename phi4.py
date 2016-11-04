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


class Phi4():
    """ main class """
    def __init__(self, m, L):

        self.m = m
        self.L = L
        self.nchunks = 10

        self.basis = {}
        self.h0 = {}
        self.V = {k:{} for k in {-1,1}}
        self.DeltaH = {}
        self.VLh = {}
        self.Vhl = {}
        self.Vll = {}
        self.V0V4 = {}
        self.V2V4 = {}
        self.V4V4 = {}
        # self.VhhDiagList = {}
        self.basisH = {}
        self.basisl = {}

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


    def genHEBasis(self, k, basisl, EL):
        """ Generate a high-energy basis from a set of tails
        k: parity quantum number
        basisl: Basis instance containing the set of tails
        EL: maximal energy of the generated basis
        """

        self.basisl[k] = basisl

        # Generate all the operators between the selected states and the states
        # in the range [0, EL]
        Vlist = V4OpsSelectedFull(basisl, EL)
        vectorset = set()

        for V in Vlist:
            for v in V.genBasis(basisl, EL):
# Don't add twice states connected by parity inversion
                if v not in vectorset and v[::-1] not in vectorset:
                    vectorset.add(v)

        helper = Vlist[0].helper

        self.basisH[k] = Basis(k, [helper.torepr1(v) for v in vectorset], helper)


    # @profile
    def computeHEVs(self, k, EL):
        """
        Compute the matrices involving the high-energy states below EL
        k: parity quantum number
        EL: local cutoff
        """

        # NOTE matrix subscript notation:
# "l": selected low-energy state
# "L": generic low-energy state
# "h": selected high-energy state

        #################################
        # Generate the Vlh matrix
        #################################
        basis = self.basisl[k]
        lookupbasis = self.basisH[k]

        Vlist = V4OpsSelectedFull(basis, EL)

        self.Vhl[k] = Matrix(basis, lookupbasis,
                self.buildMatrix(Vlist, basis, lookupbasis)*self.L)


        ##############################
        # Generate the VhL matrix
        ##############################

        print("Computing VhL...")

        basis = self.basisH[k]
        lookupbasis = self.basis[k]

        Vlist = V4OpsSelectedFull(basis, lookupbasis.Emax)

        self.VLh[k] = Matrix(basis, lookupbasis,
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

        Vlist = V6OpsSelectedHalf(basis)

        self.Vll[k][6] = Matrix(basis,basis,
                self.buildMatrix(Vlist, basis, basis, ignKeyErr=True,
                    sumTranspose=True)*self.L)

        ###################################
        # Generate all the "bilocal" matrices on the selected
        # low-energy states
        ###################################
        self.V0V4[k] = self.Vll[k][4]*self.L

        Vlist = V2V4Ops1(basis)
        self.V2V4[k] = Matrix(basis,basis,
                self.buildMatrix(Vlist,basis,basis,ignKeyErr=True,
                    sumTranspose=True, subDiag=False)*self.L**2)
        Vlist = V2V4Ops2(basis)
        self.V2V4[k] += Matrix(basis,basis,
                self.buildMatrix(Vlist,basis,basis,ignKeyErr=True,
                        sumTranspose=False)*self.L**2)

        Vlist = V4V4Ops1(basis)
        self.V4V4[k] = Matrix(basis,basis,
                self.buildMatrix(Vlist,basis,basis,ignKeyErr=True,
                    sumTranspose=True, subDiag=False)*self.L**2)
        Vlist = V4V4Ops2(basis)
        self.V4V4[k] += Matrix(basis,basis,
                self.buildMatrix(Vlist,basis,basis,ignKeyErr=True,
                        sumTranspose=False)*self.L**2)



    # @profile
    def computeDeltaH(self, k, ET, EL, eps, EL3=None, maxntails=None):
# Compute the full DeltaH = DH2 * (DH2-DH3)^-1 * DH2  matrix

        helper = self.basis[k].helper
        # Subset of the full low energy states
        subbasisL = self.basis[k].sub(lambda v: helper.energy(v)<=ET)

        VV2 = renorm.renVV2(g4=self.g4, EL=EL, eps=eps).VV2
        # Dictionary of local renormalization coefficients for the g^2 term

        if EL3==None:
            EL3 = EL

        # FIXME change this condition
        # Local + non-local
        if (EL != ET):
            basisl = self.basisl[k]
            helper = self.basisH[k].helper

            # Subset of the selected low energy states
            # print(self.basisl)
            vectorlist = [v for v in basisl if helper.energy(v)<=ET]
            if maxntails != None:
                vectorlist = vectorlist[:maxntails]
            subbasisl = Basis(k, vectorlist, basisl.helper)
            self.ntails = subbasisl.size

            # Subset of the selected high energy states
            # subbasisH = self.basisH[k].sub(lambda v: ET<helper.energy(v)<=EL)


            # Subset of the selected high energy states
            helper = self.basisH[k].helper
            vectorlist = [v for v in self.basisH[k] if ET<helper.energy(v)<=EL]
            subbasisH = Basis(k, vectorlist, helper)

            # Propagator on the high-energy states
            energyArr = array(subbasisH.energyList)
            propagator = scipy.sparse.spdiags(1/(eps-energyArr),0,subbasisH.size,
                    subbasisH.size)


            #############################
            # Construct the matrices needed to compute DeltaH
            #############################
            Vhl = self.Vhl[k].sub(subbasisl, subbasisH).M
            Vlh = Vhl.transpose()
            VLh = self.VLh[k].sub(subbasisH, subbasisL).M
            VhL = VLh.transpose()


            VlL = {}
            VLl = {}
            for n in (0,2,4):
                VlL[n] = self.V[k][n].sub(subbasisL, subbasisl).M
                VLl[n] = VlL[n].transpose()

            Vll = {}
            # for n in (0,2,4,6,8):
            for n in (0,2,4,6):
                Vll[n] = self.Vll[k][n].sub(subbasisl, subbasisl).M

            V0V4 = self.V0V4[k].sub(subbasisl,subbasisl).M
            V2V4 = self.V2V4[k].sub(subbasisl,subbasisl).M
            V4V4 = self.V4V4[k].sub(subbasisl,subbasisl).M



            #######################################
            # Construct DH2
            #######################################

            # XXX Sorry, for now the subscripts are confusing, need to sort this out
            DH2lL = VhL*propagator*Vlh*self.g4**2
            # print("DH2lL", DH2lL.M.shape)
            DH2lL += VV2[0]*VlL[0] + VV2[2]*VlL[2] + VV2[4]*VlL[4]
            DH2Ll = DH2lL.transpose()
            DH2ll = Matrix(subbasisl, subbasisL, DH2Ll).sub(subbasisl, subbasisl).M.tocsc()




###########################
# Computation of DH3
##########################

            # TODO Add the local parts
            # NOTE Trick to save memory


            # Subset of the selected high energy states
            helper = self.basisH[k].helper
            vectorlist = [v for v in self.basisH[k] if ET<helper.energy(v)<=EL3]
            subbasisH = Basis(k, vectorlist, helper)

            # Propagator on the high-energy states
            energyArr = array(subbasisH.energyList)
            propagator = scipy.sparse.spdiags(1/(eps-energyArr),0,subbasisH.size,
                    subbasisH.size)

            Vhl = self.Vhl[k].sub(subbasisl, subbasisH).M
            Vlh = Vhl.transpose()

            DH3ll = scipy.sparse.csc_matrix((subbasisl.size, subbasisl.size))

            basis = subbasisH
            Vhhlist = V4OpsSelectedHalf(basis)

            chunklen = int(math.ceil(basis.size/self.nchunks))
            idxLists = [range(basis.size)[x:x+chunklen] for x in
                    range(0, basis.size, chunklen)]

            for n in range(self.nchunks):
                ##############################
                # Generate the Vhh matrix
                ##############################
                # print("Computing Vhh chunk ", n)

                idxList = idxLists[n]

                # NOTE Trick to save memory: we never compute explicitly the full matrix Vhh
# TODO Do this more elegantly with iterators ?
                VhhHalfPart =  self.buildMatrix(Vhhlist, basis, basis,
                        ignKeyErr=True, sumTranspose=False, idxList=idxList)*self.L

                VhhDiagPart = scipy.sparse.spdiags(VhhHalfPart.diagonal(),0,basis.size,
                        basis.size)

                DH3llPart = Vhl*propagator*VhhHalfPart*propagator*Vlh*self.g4**3
                DH3llPart += DH3llPart.transpose()
                DH3llPart -= Vhl*propagator*VhhDiagPart*propagator*Vlh*self.g4**3
                DH3ll += DH3llPart

                del VhhHalfPart
                gc.collect()


            return DH2lL*scipy.sparse.linalg.inv(DH2ll-DH3ll)*DH2Ll


        # Only DH2 with local expansion
# FIXME it should not depend only on the value of EL and ET
        elif EL==ET:
            VLL = {}
            for n in (0,2,4):
                VLL[n] = self.V[k][n].sub(subbasisL, subbasisL)

            return (VV2[0]*VLL[0] + VV2[2]*VLL[2] + VV2[4]*VLL[4]).M




    def calcVV3(self, ETlist, eps):
        print("Calculating VVV renorm coefficients")
        self.VV3 =  renorm.renVV3(m=self.m, eps=eps, ETlist=ETlist)


    def setCouplings(self, g0=0, g2=0, g4=0):
        self.g0 = g0
        self.g2 = g2
        self.g4 = g4



    def computeEigval(self, k, ET, ren, EL=None, EL3=None, eps=None, neigs=10, maxntails=None):
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
        else:
            DeltaH = self.computeDeltaH(k=k, ET=ET, EL=EL, eps=eps, maxntails=maxntails, EL3=EL3)
            compH = (Hraw + DeltaH)

        self.compSize = compH.shape[0]

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


    @staticmethod
    def matrixfname(L, Emax, k, occmax):
        return "matrices/L={0:}_Emax={1:}_k={2:}_occmax={3:}".format(L,Emax,k,occmax)

    def saveMatrix(self, k, Emax):
        """ Saves the potential and free hamiltonian to file """
        L = self.helper.L
        m = self.helper.m

        fname = self.matrixfname(L, Emax, k, self.basis[k].occmax)
        t = (fname, L, m, k, \
            Emax, self.basis[k].occmax, \
            self.h0[k].M.data,self.h0[k].M.row,self.h0[k].M.col, \
            self.V[k][0].M.data,self.V[k][0].M.row,self.V[k][0].M.col, \
            self.V[k][2].M.data,self.V[k][2].M.row,self.V[k][2].M.col, \
            self.V[k][4].M.data,self.V[k][4].M.row,self.V[k][4].M.col, \
            )
        scipy.savez(*t)

    def loadMatrix(self, L, Emax, k, occmax):
        """ Loads the potential and free hamiltonian from file """
        fname = self.matrixfname(L, Emax, k, occmax)+".npz"
        f = scipy.load(fname)
        L = f['arr_0'].item()
        m = f['arr_1'].item()

        k = f['arr_2'].item()
        Emax = f['arr_3'].item()
        occmax = f['arr_4'].item()

        self.buildBasis(L=L, Emax=Emax, m=m, k=k, occmax=occmax)
        basis = self.basis[k]

        z = 5
        self.h0[k] = Matrix(basis, basis, scipy.sparse.coo_matrix(
            (f['arr_'+str(z)], (f['arr_'+str(z+1)],f['arr_'+str(z+2)])),
            shape=(basis.size, basis.size)))
        z = 8
        for n in (0,1,2):
            self.V[k][2*n] = Matrix(basis, basis, scipy.sparse.coo_matrix(
            (f['arr_'+str(z+3*n)],(f['arr_'+str(z+1+3*n)], f['arr_'+str(z+2+3*n)])),
              shape=(basis.size, basis.size)))
