import scipy
import scipy.sparse.linalg
import scipy.sparse
import math
from math import factorial
from statefuncs import Basis, helper
from oscillators import V2Operators, V4Operators, V4OperatorsLH
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
        self.DeltaH2 = {}
        self.VLH = {}
        self.basisH = {}

        self.eigenvalues = {"raw":{}, "renlocal":{}}
        self.eigenvectors = {"raw":{}, "renlocal":{}}

        self.compBasisSize = {1: None, -1:None}
        # Holds basis sizes used to actually compute the eigenvalues

        scipy.set_printoptions(precision=15)


    def buildBasis(self, Emax, occmax=None):
        """ Builds the full Hilbert space basis """
        self.helper = helper(self.m, self.L, Emax)

        self.basis = Basis.fromScratch(Emax=Emax, helper=self.helper, occmax=occmax)


    def buildMatrix(self):
        """ Builds the full hamiltonian in the basis of the free hamiltonian.
        This is computationally intensive.
        It can be skipped by loading the matrix from file """


        Vlist = {2:V2Operators(self.helper, self.basis),
                    4:V4Operators(self.helper, self.basis)}
        L = self.helper.L


        for k in (-1,1):
            basis = self.basis[k]
            lookupbasis = self.basis[k]
            Emax = basis.Emax

            self.h0[k] = Matrix(basis, basis,
                    scipy.sparse.spdiags([self.helper.energy(v) for v in basis],
                        0,basis.size,basis.size)).to('coo')

            # Construct the quadratic and quartic potential matrices
            for n in (2,4):
                # Will construct the sparse matrix in the COO format
                data = []
                row = []
                col = []

                for V in Vlist[n]:
                    computeME = V.computeMatrixElements
                    for i in range(basis.size):
                        colpart, datapart = computeME(basis,i,lookupbasis)
                        data += datapart
                        col += colpart
                        row += [i]*len(colpart)


                V = scipy.sparse.coo_matrix((data,(row,col)), shape=(basis.size,basis.size))

                diag_V = scipy.sparse.spdiags(V.diagonal(),0,basis.size,basis.size)
                # This should resum elements in the same coordinate
                self.V[k][n] = Matrix(basis,basis,V+V.transpose()-diag_V).to('coo')*L


            # Construct the identity potential matrix
            idM = scipy.sparse.eye(basis.size)
            self.V[k][0] = Matrix(basis, basis, idM).to('coo')*L

    def computeDH2(self, k, subbasis, Emin, ET):

        params = helper(self.helper.m, self.helper.L,ET)
        L = self.helper.L

        # Generate the high-energy basis
        Vlist = V4OperatorsLH(params, subbasis, Emin, ET)
        vectorset = set()


        for V in Vlist:
            for v in V.genBasis(subbasis, Emin, ET):
                if v not in vectorset and v[::-1] not in vectorset:
                    vectorset.add(v)

        self.basisH[k] = Basis(k, [params.torepr1(v) for v in vectorset], params)


        # Generate the matrix
        basis = subbasis
        lookupbasis = self.basisH[k]

        data = []
        row = []
        col = []

        for V in Vlist:
            computeME = V.computeMatrixElements
            for i in range(basis.size):
                colpart, datapart = computeME(basis,i,lookupbasis)
                data += datapart
                col += colpart
                row += [i]*len(colpart)

        V = scipy.sparse.coo_matrix((data,(row,col)), shape=(basis.size,lookupbasis.size))
        self.VLH[k] = Matrix(basis,lookupbasis,V).to('coo')*L


    def setCouplings(self, g0, g2, g4):
        self.g0 = g0
        self.g2 = g2
        self.g4 = g4
        c = 2.*sum([1/(2.*pi)*scipy.special.kn(0,n*m*L) for n in range(1,10)])
        self.m1 = m*exp(-2.*pi*c)


    def renlocal(self, Emax, Er):
        self.g0r, self.g2r, self.g4r = \
            renorm.renlocal(self.g0, self.g2, self.g4, Emax, m=self.m1, Er=Er)


    def computeHamiltonian(self, k, Emax, ren):
        if ren=="raw":
            V = self.V[k][0]*self.g0 + self.V[k][2]*self.g2 + self.V[k][4]*self.g4
        elif ren=="renlocal":
            V = self.V[k][0]*self.g0r + self.V[k][2]*self.g2r + self.V[k][4]*self.g4r
        else:
            raise ValueError()

        H0 = self.h0[k]
        H = H0 + V

        self.compH = H

        # basisL = Basis.fromBasis(self.basis[k], lambda v: v.energy <= Emax)
        # basisH = Basis.fromBasis(self.basis[k], lambda v: v.energy > Emax)
        # Hll = H.sub(basisL, basisL).M
        # gramL = scipy.sparse.eye(basisL.size)

        # self.compH = Hll
        # self.gram = gramL

    def computeEigval(self, ren, k, sigma=0, neigs=10):
        """ Sets the internal variables self.eigenvalues
        k : field parity
        Emax : max energy of truncated Hilbert space
        ren : renormalization procedure "raw" or "renlocal"
        neigs : number of eigenvalues to compute
        sigma : point around which we should look for eigenvalues.
        """

        (self.eigenvalues[ren][k], eigenvectorstranspose) = \
            scipy.sparse.linalg.eigsh(self.compH, M=self.gram, k=neigs, sigma=sigma,
                            which='LM', return_eigenvectors=True)
        self.eigenvectors[ren][k] = eigenvectorstranspose.T

        gc.collect()

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
        else:
            raise ValueError("ren value not valid")


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
