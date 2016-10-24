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
from matrix import Matrix
from scipy import exp, pi, array
from pathos.multiprocessing import ProcessingPool as Pool
from multiprocessing import Process, Queue


class Phi4():
    """ main class """
    def __init__(self, m, L):

        self.m = m
        self.L = L
        self.nproc = 3

        self.basis = {}
        self.h0 = {}
        self.V = {k:{} for k in {-1,1}}
        self.DeltaH = {}
        self.VLh = {}
        self.Vhl = {}
        self.VhhHalfList = {}
        self.VhhDiagList = {}
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


    def buildMatrixChunk(self, Vlist, basis, lookupbasis, idxList, statePos,
            helper, ignKeyErr=False):

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

        # This should resum duplicate entries
        return V.tocsc()


    # @profile
    def buildMatrix(self, Vlist, basis, lookupbasis, ignKeyErr=False, sumTranspose=False,
            parallel=False):
        """
        Compute a potential matrix from the operators in Vlist between two bases
        Vlist: list of Operator instances (for instance the
        a^a^a^a^, a^a^a^a and a^a^aa and parts of V4)
        basis: set of states corresponding to the row indexes of the matrix
        lookupbasis: set of states corresponding to column indexes of the matrix
        ignKeyErr: this should be set to True only for Vhh (because some states generated
        don't belong to the basis)
        """

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



        if parallel==False:
            idxList = range(basis.size)

            # Build a single sparse matrix in the CSC format
            V = self.buildMatrixChunk(Vlist, basis, lookupbasis, idxList, statePos, helper,
                    ignKeyErr=ignKeyErr)

            if sumTranspose:
                # Add the matrix to its transpose and subtract the diagonal
                diag_V = scipy.sparse.spdiags(V.diagonal(),0,basis.size,basis.size).tocsc()
                return (V+V.transpose()-diag_V)*self.L
            else:
                return V*self.L

        else:
            # Split the basis index list into chunks, and generate a list of separate
            # sparse matrices in the CSC format. Not adding them up saves memory.
            chunklen = int(math.ceil(basis.size/self.nproc))
            idxLists = [range(basis.size)[x:x+chunklen] for x in
                    range(0, basis.size, chunklen)]

            def f(q, idxList):
                q.put(self.buildMatrixChunk(Vlist, basis, lookupbasis,
                            idxList, statePos, helper, ignKeyErr=ignKeyErr)*self.L)


            q = Queue()
            p = {}
            V = []
            for n in range(self.nproc): p[n] = Process(target=f, args=(q,idxLists[n]))
            for n in range(self.nproc): p[n].start()
            for n in range(self.nproc): V.append(q.get())
            return V



            # with Pool(self.nproc) as p:
                # return p.map(f, idxLists)


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
                    self.buildMatrix(Vlist[n], basis, basis, sumTranspose=True))

        # Construct the identity potential matrix
        idM = scipy.sparse.eye(basis.size)
        self.V[k][0] = Matrix(basis, basis, idM)


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
    def computeHEVs(self, k, EL, EL3=None):
        """
        Compute the matrices involving the high-energy states below EL
        k: parity quantum number
        EL: local cutoff
        """

        # NOTE matrix subscript notation:
# "l": selected low-energy state
# "L": generic low-energy state
# "h": selected high-energy state

        # We might want to choose EL3 < EL because Vhh is expensive to generate
        if EL3 == None:
            EL3 = EL


        #################################
        # Generate the Vlh matrix
        #################################
        basis = self.basisl[k]
        lookupbasis = self.basisH[k]

        Vlist = V4OpsSelectedFull(basis, EL)

        self.Vhl[k] = Matrix(basis, lookupbasis,
                self.buildMatrix(Vlist, basis, lookupbasis))


        ##############################
        # Generate the VhL matrix
        ##############################


        print("Computing VhL...")

        basis = self.basisH[k]
        lookupbasis = self.basis[k]

        Vlist = V4OpsSelectedFull(basis, lookupbasis.Emax)

        self.VLh[k] = Matrix(basis, lookupbasis,
                self.buildMatrix(Vlist, basis, lookupbasis))


        ##############################
        # Generate the Vhh matrix
        ##############################

        print("Computing Vhh...")

        basis = self.basisH[k]

        Vlist = V4OpsSelectedHalf(basis)

        # NOTE Trick to save memory: we never compute explicitly the full matrix Vhh
        self.VhhHalfList[k] =  self.buildMatrix(Vlist, basis, basis, ignKeyErr=True,
                    sumTranspose=False, parallel=True)

        self.VhhDiagList[k] = \
            [scipy.sparse.spdiags(VhhHalfChunk.diagonal(),0,basis.size,basis.size)
                for VhhHalfChunk in self.VhhHalfList[k]]


    # @profile
    def computeDeltaH(self, k, ET, EL, eps, maxntails):
# Compute the full DeltaH = DH2 * (DH2-DH3)^-1 * DH2  matrix

        helper = self.basis[k].helper
        # Subset of the full low energy states
        subbasisL = self.basis[k].sub(lambda v: helper.energy(v)<=ET)

        # Dictionary of local renormalization coefficients
        VV2 = renorm.renlocal(g4=self.g4, EL=EL, m=self.m, eps=eps).VV2

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


            # Propagator and projector on the high-energy states
            basisH = self.basisH[k]
            propVec = []
            for e in basisH.energyList:
                if ET < e <= EL:
                    propVec.append(1/(eps-e))
                else:
                    propVec.append(0)
            propagator = scipy.sparse.spdiags(propVec,0,basisH.size,basisH.size)


            #############################
            # Construct DH
            #############################
            Vhl = self.Vhl[k].sub(subbasisl, basisH).M
            Vlh = Vhl.transpose()
            VLh = self.VLh[k].sub(basisH, subbasisL).M
            VhL = VLh.transpose()
            # Vhh = self.Vhh[k].sub(subbasisH, subbasisH)
            VhhHalfList = self.VhhHalfList[k]
            VhhDiagList = self.VhhDiagList[k]


            VlL = {}
            VLl = {}
            for n in (0,2,4):
                VlL[n] = self.V[k][n].sub(subbasisL, subbasisl).M
                VLl[n] = VlL[n].transpose()

            Vll = {}
            # for n in (0,2,4,6,8):
            for n in (0,2,4):
                Vll[n] = self.V[k][n].sub(subbasisl, subbasisl).M


            # XXX Sorry, for now the subscripts are confusing, need to sort this out
            DH2lL = VhL*propagator*Vlh*self.g4**2
            # print("DH2lL", DH2lL.M.shape)
            DH2lL += VV2[0]*VlL[0] + VV2[2]*VlL[2] + VV2[4]*VlL[4]
            DH2Ll = DH2lL.transpose()
            DH2ll = Matrix(subbasisL, subbasisl, DH2Ll).sub(subbasisl, subbasisl).M.tocsc()

            # TODO Add the local parts
            # NOTE Trick to save memory
            DH3ll = scipy.sparse.csc_matrix((subbasisl.size, subbasisl.size))
            for VhhHalfChunk, VhhDiagChunk in zip(VhhHalfList, VhhDiagList):
                DH3llPart = Vhl*propagator*VhhHalfChunk*propagator*Vlh*self.g4**3
                DH3llPart += DH3llPart.transpose()
                DH3llPart -= Vhl*propagator*VhhDiagChunk*propagator*Vlh*self.g4**3
                DH3ll += DH3llPart

            return DH2lL*scipy.sparse.linalg.inv(DH2ll-DH3ll)*DH2Ll


        # Only local
        elif EL==ET:
            VLL = {}
            for n in (0,2,4):
                VLL[n] = self.V[k][n].sub(subbasisL, subbasisL)

            return VV2[0]*VLL[0] + VV2[2]*VLL[2] + VV2[4]*VLL[4]




    def setCouplings(self, g0=0, g2=0, g4=0):
        self.g0 = g0
        self.g2 = g2
        self.g4 = g4



    def computeEigval(self, k, ET, ren, EL=None, eps=None, neigs=10, maxntails=None):
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

        Hraw = (self.h0[k] + self.g4*self.V[k][4]).sub(subbasisL, subbasisL)

        if ren=="raw":
            compH = Hraw.M
        else:
            DeltaH = self.computeDeltaH(k, ET, EL, eps, maxntails)
            compH = (Hraw + DeltaH).M

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
        eigs = self.eigenvalues[ren]
        # Subtract vacuum energies
        if k==1:
            return scipy.array([x-self.vacuumE(ren=ren) for x in eigs[k][1:]])
        elif k==-1:
            return scipy.array([x-self.vacuumE(ren=ren) for x in eigs[k]])


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
