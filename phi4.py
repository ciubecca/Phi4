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



class Phi4():
    """ main FVHT class """
    def __init__(self, m, L):

        self.m = m
        self.L = L

        self.basis = {}
        self.h0 = {}
        self.V = {k:{} for k in {-1,1}}
        self.DeltaH = {}
        self.VLh = {}
        self.Vhl = {}
        self.Vhh = {}
        self.basisH = {}
        self.basisl = {}

        self.eigenvalues = {"raw":{}, "ren":{}}
        self.eigenvectors = {"raw":{}, "ren":{}}

        self.compSize = {}
        # Holds basis sizes used to actually compute the eigenvalues

        scipy.set_printoptions(precision=15)


    def buildBasis(self, Emax, occmax=None):
        """ Builds the full Hilbert space basis """
        self.basis = Basis.fromScratch(m=self.m, L=self.L, Emax=Emax, occmax=occmax)


    def buildMatrix(self, Vlist, basis, lookupbasis, ignKeyErr=False):
        """
        Compute a potential matrix from the operators in Vlist
        between two bases
        """

        if basis.helper.nmax > lookupbasis.helper.nmax:
            helper = basis.helper
        else:
            helper = lookupbasis.helper

        # Contains also the P-reversed states
        # NOTE: using arrays is much less efficient!
        statePos = {}
        for i,state in enumerate(lookupbasis.stateList):
            statePos[tuple(helper.torepr2(state))] = i
            statePos[tuple(helper.torepr2(state)[::-1])] = i

        # Will construct the sparse matrix in the COO format
        data = []
        row = []
        col = []

        for V in Vlist:
            for i in range(basis.size):
                colpart, datapart = \
                    V.computeMatrixElements(basis,i,lookupbasis,
                            statePos=statePos, helper=helper, ignKeyErr=ignKeyErr)
                data += datapart
                col += colpart
                row += [i]*len(colpart)

        V = scipy.sparse.coo_matrix((data,(row,col)),
                shape=(basis.size,lookupbasis.size))

        if basis==lookupbasis:
        # If the two bases are equal, assumes that only half of the matrix was computed
            diag_V = scipy.sparse.spdiags(V.diagonal(),0,basis.size,basis.size)
            # This should resum duplicate entries
            return (V+V.transpose()-diag_V).tocsc()
        else:
            # XXX This should resum duplicate entries
            return V.tocsc()



    def computePotential(self, k):
        """
        Builds the potential matrices and the free Hamiltonian.
        """

        basis = self.basis[k]

        self.h0[k] = Matrix(basis, basis,
                scipy.sparse.spdiags(basis.energyList,0,basis.size,basis.size))

        Vlist = {2:V2OpsHalf(basis), 4:V4OpsHalf(basis)}
        for n in (2,4):
            self.V[k][n] = Matrix(basis,basis,
                    self.buildMatrix(Vlist[n], basis, basis))*self.L

        # Construct the identity potential matrix
        idM = scipy.sparse.eye(basis.size)
        self.V[k][0] = Matrix(basis, basis, idM)*self.L

    # @profile
    def genHEBasis(self, k, basisl, EL):

        self.basisl[k] = basisl

        # Generate all the operators between the selected states and the states
        # between in the range [0, EL]
        Vlist = V4OpsSelectedFull(basisl, EL)
        vectorset = set()

        for V in Vlist:
            for v in V.genBasis(basisl, EL):
                if v not in vectorset and v[::-1] not in vectorset:
                    vectorset.add(v)

        helper = Vlist[0].helper

        self.basisH[k] = Basis(k, [helper.torepr1(v) for v in vectorset], helper)

    # @profile
    def computeHEVs(self, k, EL, EL3=None):
        # NOTE "l" denotes a selected low-energy state, while "L" a
        # generic low-energy state

        # We might want to choose EL3 < EL because Vhh is expensive to generate
        if EL3 == None:
            EL3 = EL


        # NOTE Notation:
        # SL : Selected Low
        # SH : Selected High
        # FL : Full Low

        #################################
        # Generate the SL-SH matrix
        #################################
        basis = self.basisl[k]
        lookupbasis = self.basisH[k]

        Vlist = V4OpsSelectedFull(basis, EL)

        self.Vhl[k] = Matrix(basis, lookupbasis,
                self.buildMatrix(Vlist, basis, lookupbasis)*self.L)


        ##############################
        # Generate the SH-FL matrix
        ##############################
        basis = self.basisH[k]
        lookupbasis = self.basis[k]

        Vlist = V4OpsSelectedFull(basis, lookupbasis.Emax)

        self.VLh[k] = Matrix(basis, lookupbasis,
                self.buildMatrix(Vlist, basis, lookupbasis)*self.L)


        ##############################
        # Generate the SH-SH matrix
        ##############################
        basis = self.basisH[k]

        Vlist = V4OpsSelectedHalf(basis)

        self.Vhh[k] = Matrix(basis, basis,
                self.buildMatrix(Vlist, basis, basis, ignKeyErr=True)*self.L)


    def computeDeltaH(self, k, ET, EL, eps):

        helper = self.basisH[k].helper

        # Subset of the full low energy states
        subbasisL = self.basis[k].sub(lambda v: helper.energy(v)<=ET)

        # Subset of the selected low energy states
        subbasisl = self.basisl[k].sub(lambda v: helper.energy(v)<=ET)
        self.ntails = subbasisl.size

        # Subset of the selected high energy states
        subbasisH = self.basisH[k].sub(lambda v: ET<helper.energy(v)<=EL)

        #############################
        # Construct DH
        #############################
        Vhl = self.Vhl[k].sub(subbasisl, subbasisH)
        Vlh = Vhl.transpose()
        VLh = self.VLh[k].sub(subbasisH, subbasisL)
        VhL = VLh.transpose()
        Vhh = self.Vhh[k].sub(subbasisH, subbasisH)

        propagator = Matrix(subbasisH, subbasisH,
                scipy.sparse.spdiags([1/(eps-e) for e in subbasisH.energyList],
                0, subbasisH.size, subbasisH.size))

        # Dictionary of local renormalization coefficients
        deltag = renorm.renlocal(g4=self.g4, EL=EL, m=self.m, eps=eps)

        VLl = {}
        for n in (0,2,4):
            VLl[n] = self.V[k][n].sub(subbasisl, subbasisL)

        Vll = {}
        # for n in (0,2,4,6,8):
        for n in (0,2,4):
            Vll[n] = self.V[k][n].sub(subbasisl, subbasisl)


        # The first index in deltag denotes the order in V
        # The second index is a tuple denoting the type of local operators
        DH2Ll = Vhl*propagator*VLh*self.g4**2 \
                + deltag[2][0]*VLl[0]+deltag[2][2]*VLl[2]+deltag[2][4]*VLl[4]
        DH2lL = DH2Ll.transpose()
        DH2ll = DH2Ll.sub(subbasisl, subbasisl)

        # The first index in deltag denotes the order in V
        # The second index is a tuple denoting the type of local operators
        # TODO Add bi-local operators. Do these operators need to be generated separately?
        # DH3ll = Vhl*propagator*Vhh*propagator*Vlh*self.g4**3 +
                # + deltag[3][0]*Vll[0]+deltag[3][2]*Vll[2]+deltag[3][4]*Vll[4]
                # + deltag[3][6]*Vll[6]+deltag[3][8]*Vll[8]

        # TODO Add the local parts
        DH3ll = Vhl*propagator*Vhh*propagator*Vlh*self.g4**3

        # print(DH2lL.shape, DH2ll.shape, DH2Ll.shape)
        # print(scipy.sparse.linalg.inv(DH2ll).shape)

        return DH2lL*(DH2ll-DH3ll).inverse()*DH2Ll




    def setCouplings(self, g0, g2, g4):
        self.g0 = g0
        self.g2 = g2
        self.g4 = g4
        # c = 2.*sum([1/(2.*pi)*scipy.special.kn(0,n*m*L) for n in range(1,10)])
        # self.m1 = m*exp(-2.*pi*c)



    def computeEigval(self, k, ET, ren, EL=None, eps=None, neigs=10):

        # Subset of low energy states
        helper = self.basis[k].helper

        subbasisL = self.basis[k].sub(lambda v: helper.energy(v)<=ET)
        self.subbasisL = subbasisL

        Hraw = (self.h0[k] + self.g4*self.V[k][4]).sub(subbasisL, subbasisL)

        if ren=="raw":
            compH = Hraw.M
        elif ren=="ren":
            DeltaH = self.computeDeltaH(k, ET, EL, eps)
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
