import sys
import scipy
import scipy.sparse.linalg
import scipy.sparse
import scipy.interpolate
import math
from operator import attrgetter
import gc
import copy
import statefuncs
from math import factorial
from statefuncs import Basis, NotInBasis, omega, State
import oscillators
from oscillators import NormalOrderedOperator as NOO
import collections
import renorm
import itertools
import finiteVolH
from scipy import exp, pi


# NOT IMPORTANT. I use this "tolerance" parameter throughout the code to avoid possible numerical issues when confronting energies with Emax
tol = 0.0001

""" P denotes spatial parity, while K field parity. For now only the P-even sector is implemented """

def comb(*x):
    """ computes combinatorial factor for list of elements """
    return factorial(len(x))/scipy.prod(list(map(factorial,collections.Counter(x).values())))

class Matrix():
    """ Matrix with specified state bases for row and column indexes.
    This class is useful to easily extract submatrices """
    def __init__(self, basisI, basisJ, M=None):
        self.basisI = basisI
        self.basisJ = basisJ

        if(M == None):
            self.M = scipy.sparse.coo_matrix((basisI.size, 1))
        else:
            self.M = M
            self.check()

    def addColumn(self, newcolumn):
        m = scipy.sparse.coo_matrix(newcolumn).transpose()
        self.M = scipy.sparse.hstack([self.M,m])

    def finalize(self):
        self.M = self.M.tocsc()[:,1:].tocoo()
        self.check()

    def check(self):
        if self.M.shape != (self.basisI.size, self.basisJ.size):
            raise ValueError('Matrix shape inconsistent with given bases')

    def __add__(self, other):
        """ Sum of matrices """

        return Matrix(self.basisI, self.basisJ, self.M+other.M)

    def __mul__(self, other):
        """ Multiplication of matrix with matrix or number"""
        if(other.__class__ == self.__class__):
            return Matrix(self.basisI, other.basisJ, self.M*other.M)
        else:
            return Matrix(self.basisI, self.basisJ, self.M*float(other))

    def to(self, form):
        """ Format conversion """
        return Matrix(self.basisI, self.basisJ, self.M.asformat(form))

    def sub(self, subBasisI=None, subBasisJ=None):
        """ This extracts a submatrix given a subspace of the initial vector space, both for rows and columns """

        if subBasisI != None and subBasisJ != None:
            rows = [self.basisI.lookup(state)[1]  for state in subBasisI]
            columns = [self.basisJ.lookup(state)[1]  for state in subBasisJ]
            return Matrix(subBasisI, subBasisJ, self.M.tocsr()[scipy.array(rows)[:,scipy.newaxis],scipy.array(columns)])

        elif subBasisI != None and subBasisJ == None:
            rows = [self.basisI.lookup(state)[1]  for state in subBasisI]
            return Matrix(subBasisI, self.basisJ, self.M.tocsr()[scipy.array(rows),:])

        elif subBasisI == None and subBasisJ != None:
            columns = [self.basisJ.lookup(state)[1]  for state in subBasisJ]
            return Matrix(self.basisI, subBasisJ, self.M.tocsr()[:,scipy.array(columns)])

        else:
            return self

    def transpose(self):
        return Matrix(self.basisJ, self.basisI, self.M.transpose())

class Phi1234():
    """ main FVHT class """
    def __init__(self):
        self.L=None
        self.m=None

        self.h0 = {1 : None, -1 : None}
        self.potential = { 1 : { }, -1 : { }}
        self.h0Sub = {1 : None, -1 : None}
        self.H = {1 : None, -1 : None}
        self.V = {1 : { }, -1 : { }}

        self.eigenvalues = {1: None, -1:None}
        self.eigsrenlocal = {1: None, -1:None}
        self.eigsrensubl = {1: None, -1:None}
        self.eigenvectors = {1: None, -1:None}
        # Eigenvalues and eigenvectors for different K-parities

        self.compBasisSize = {1: None, -1:None}
        # Holds basis sizes used to actually compute the eigenvalues

        scipy.set_printoptions(precision=15)

        self.basis = {1: None, -1:None}
        self.fullBasis = {1: None, -1:None}

    def buildFullBasis(self,k,L,m,Emax):
        """ Builds the full Hilbert space basis """

        self.L = float(L)
        self.m = float(m)

        self.fullBasis[k] = Basis(L=self.L, Emax=Emax, m=self.m, K=k)


    def buildBasis(self, k, Emax):
        """ Builds the Hilbert space basis for which the Hamiltonian to actually diagonalize
        is calculated (in general it's a subspace of fullBasis) """

        self.basis[k] = Basis(m=self.m, L=self.L, Emax=Emax, k=k, nmax=self.fullBasis[k].nmax)
        # We use the vector length (nmax) of the full basis. In this way we can compare elements between the two bases
        # print('nmax :', self.basis[k].nmax)
        self.Emax = float(Emax)

        for nn in (0,2,4):
            self.V[k][nn] = self.potential[k][nn].sub(self.basis[k], self.basis[k]).M.tocoo()

        self.h0Sub[k] = self.h0[k].sub(self.basis[k],self.basis[k]).M.tocoo()

    def buildMatrix(self):
        """ Builds the full hamiltonian in the basis of the free hamiltonian.
        This is computationally intensive. It can be skipped by loading the matrix from file """
        L=self.L
        m=self.m

        for k in (1,-1):
            basis = self.fullBasis[k]
            lookupBasis = self.fullBasis[k]
            Emax = basis.Emax
            nmax = basis.nmax

            diagOps = {0: None, 2:None, 4:None}
            offdiagOps = {0: None, 2:None, 4:None}

            diagOps[0] = [ NOO([],[],L,m) ]

            offdiagOps[0] = []

            diagOps[2] = [ NOO([a],[a],L,m, extracoeff = 2.) for a in range(-nmax,nmax+1) ]

            offdiagOps[2] = [ NOO([a,-a],[],L,m,extracoeff=comb(a,-a))
                    for a in range(-nmax,nmax+1) if a<=-a<=nmax and
                    omega(a,L,m)+omega(-a,L,m) <= Emax+tol]

            diagOps[4] = [ NOO([a,b],[c,a+b-c],L,m, extracoeff = 6. * comb(a,b)*comb(c,a+b-c))
                    for a in range(-nmax,nmax+1) for b in range (a,nmax+1)
                    for c in range(-nmax,nmax+1) if
                    ( c<=a+b-c<=nmax
                    and (a,b) == (c,a+b-c)
                    and -Emax-tol <= omega(a,L,m)+omega(b,L,m) - omega(c,L,m)-omega(a+b-c,L,m) <=Emax+tol)]

            offdiagOps[4] = [ NOO([a,b,c,-a-b-c],[],L,m,extracoeff=comb(a,b,c,-a-b-c))
                    for a in range(-nmax,nmax+1) for b in range (a,nmax+1)
                    for c in range(b,nmax+1) if c<=-a-b-c<=nmax and
                    omega(a,L,m)+omega(b,L,m) + omega(c,L,m)+omega(-a-b-c,L,m)<= Emax+tol]  \
                + [ NOO([a,b,c],[a+b+c],L,m, extracoeff = 4. * comb(a,b,c))
                    for a in range(-nmax, nmax+1) for b in range (a,nmax+1)
                    for c in range(b,nmax+1) if
                    (-nmax<=a+b+c<=nmax
                    and -Emax-tol <= omega(a,L,m)+omega(b,L,m)+ omega(c,L,m)-omega(a+b+c,L,m) <=Emax+tol)] \
                + [ NOO([a,b],[c,a+b-c],L,m, extracoeff = 6. * comb(a,b)*comb(c,a+b-c))
                    for a in range(-nmax,nmax+1) for b in range (a,nmax+1)
                    for c in range(-nmax,nmax+1) if
                    ( c<=a+b-c<=nmax
                    and (a,b) != (c,a+b-c)
                    and sorted([abs(a),abs(b)]) < sorted([abs(c),abs(a+b-c)])
                    and -Emax-tol <= omega(a,L,m)+omega(b,L,m)- omega(c,L,m)-omega(a+b-c,L,m) <=Emax+tol)]


            #print "Number of operators:", sum([len(x) for x in offdiagOps.values()]+[len(x) for x in diagOps.values()])

            self.h0[k] = Matrix(lookupBasis, basis)
            for j in range(basis.size):
                newcolumn = scipy.zeros(lookupBasis.size)
                newcolumn[j] = basis[j].energy
                self.h0[k].addColumn(newcolumn)
            self.h0[k].finalize()

            for n in offdiagOps.keys():

                offdiag_V = Matrix(lookupBasis, basis)
                diagonal = scipy.zeros(basis.size)

                for j in range(basis.size):

                    newcolumn = scipy.zeros(lookupBasis.size)
                    for op in offdiagOps[n]:
                        try:

                            (x,i) = op.apply(basis,j,lookupBasis)

                            if(i != None):
                                newcolumn[i]+=x
                        except NotInBasis:
                            pass

                    offdiag_V.addColumn(newcolumn)

                    for op in diagOps[n]:

                        (x,i) = op.apply(basis,j,lookupBasis)
                        # It should be j=i

                        if i!= None:
                            if i != j:
                                raise RuntimeError('Non-diagonal operator')

                            diagonal[i]+=x

                offdiag_V.finalize()
                diag_V = scipy.sparse.spdiags(diagonal,0,basis.size,basis.size)

                self.potential[k][n] = (offdiag_V+offdiag_V.transpose()+Matrix(lookupBasis, basis, diag_V)).to('coo')*self.L


    def saveMatrix(self, fname):
        """ Saves the potential and free hamiltonian to file """
        t = (fname, self.L, self.m, \
            self.fullBasis[1].Emax, self.fullBasis[1].nmax,  \
            self.fullBasis[-1].Emax, self.fullBasis[-1].nmax, \
            self.h0[1].M.data,self.h0[1].M.row,self.h0[1].M.col, \
            self.potential[1][0].M.data,self.potential[1][0].M.row,self.potential[1][0].M.col, \
            self.potential[1][2].M.data,self.potential[1][2].M.row,self.potential[1][2].M.col, \
            self.potential[1][4].M.data,self.potential[1][4].M.row,self.potential[1][4].M.col, \
            self.h0[-1].M.data,self.h0[-1].M.row,self.h0[-1].M.col, \
            self.potential[-1][0].M.data,self.potential[-1][0].M.row,self.potential[-1][0].M.col, \
            self.potential[-1][2].M.data,self.potential[-1][2].M.row,self.potential[-1][2].M.col, \
            self.potential[-1][4].M.data,self.potential[-1][4].M.row,self.potential[-1][4].M.col \
            )
        scipy.savez(*t)

    def loadMatrix(self, fname):
        """ Loads the potential and free hamiltonian from file """
        f = scipy.load(fname)
        self.L = f['arr_0'].item()
        self.m = f['arr_1'].item()

        Emax = {1:f['arr_2'].item(), -1:f['arr_4'].item()}
        nmax = {1:f['arr_3'].item(), -1:f['arr_5'].item()}

        for i, k in enumerate((1,-1)):
            n = 12
            z = 6

            self.buildFullBasis(L=self.L, Emax=Emax[k], m=self.m, k=k)

            basisI = self.fullBasis[k]
            basisJ = self.fullBasis[k]

            self.h0[k] = Matrix(basisI, basisJ, scipy.sparse.coo_matrix((f['arr_'+(str(z+i*n))], (f['arr_'+(str(z+1+i*n))], f['arr_'+(str(z+2+i*n))])), shape=(basisI.size, basisJ.size)))
            self.potential[k][0] = Matrix(basisI, basisJ, scipy.sparse.coo_matrix((f['arr_'+(str(z+3+i*n))], (f['arr_'+(str(z+4+i*n))], f['arr_'+(str(z+5+i*n))])), shape=(basisI.size, basisJ.size)))
            self.potential[k][2] = Matrix(basisI, basisJ, scipy.sparse.coo_matrix((f['arr_'+(str(z+6+i*n))], (f['arr_'+(str(z+7+i*n))], f['arr_'+(str(z+8+i*n))])), shape=(basisI.size, basisJ.size)))
            self.potential[k][4] = Matrix(basisI, basisJ, scipy.sparse.coo_matrix((f['arr_'+(str(z+9+i*n))], (f['arr_'+(str(z+10+i*n))], f['arr_'+(str(z+11+i*n))])), shape=(basisI.size, basisJ.size)))

    def setCouplings(self, g0, g2, g4):
        self.g0 = g0
        self.g2 = g2
        self.g4 = g4
        c = 2.*sum([1/(2.*pi)*scipy.special.kn(0,n*self.m*self.L) for n in range(1,10)])
        self.m1 = self.m*exp(-2.*pi*c)

    def renlocal(self,Er):
        self.g0r, self.g2r, self.g4r = renorm.renlocal(self.g0,self.g2,self.g4,self.Emax,m=self.m1,Er=Er)
        self.Er = Er

    def computeHamiltonian(self, k=1, ren = False):
        """ Computes the (renormalized) Hamiltonian to diagonalize
        k : K-parity quantum number
        """

        if not(ren):
            self.H[k] = self.h0Sub[k] + self.V[k][0]*self.g0 + self.V[k][2]*self.g2 + self.V[k][4]*self.g4
        else:
            self.H[k] = self.h0Sub[k] + self.V[k][0]*self.g0r + self.V[k][2]*self.g2r + self.V[k][4]*self.g4r

        self.compBasisSize[k]=self.H[k].shape[0]


    def computeEigval(self, k=1, sigma=0, n=10, ren=False, corr=False, cutoff=None):
        """ Sets the internal variables self.eigenvalues
        k : K-parity quantum number
        n : number of eigenvalues to compute
        sigma : point around which we should look for the eigenvalues."""

        if not ren:
            (self.eigenvalues[k], eigenvectorstranspose) = scipy.sparse.linalg.eigsh(self.H[k], k=n, sigma=sigma,
                            which='LM', return_eigenvectors=True)
        else:
            (self.eigsrenlocal[k], eigenvectorstranspose) = scipy.sparse.linalg.eigsh(self.H[k], k=n, sigma=sigma,
                            which='LM', return_eigenvectors=True)
        self.eigenvectors[k] = eigenvectorstranspose.T

        if corr:
            #print "Adding subleading corrections to k=",k, " eigenvalues"

            self.eigsrensubl[k] = scipy.zeros(n)
            self.cutoff = cutoff

            for i in range(n):
                cbar = self.eigenvectors[k][i]
                if abs(sum([x*x for x in cbar])-1.0) > 10**(-13):
                    raise RuntimeError('Eigenvector not normalized')

                Ebar = self.eigsrenlocal[k][i]
                self.eigsrensubl[k][i] += Ebar
                #print self.g2, self.g4, Ebar, self.Emax, self.Er

                ktab, rentab = renorm.rensubl(self.g2, self.g4, Ebar, self.Emax, self.Er, m=self.m1, cutoff=cutoff*self.m)

                tckren = { }
                tckren[0] = scipy.interpolate.interp1d(ktab,rentab.T[0],kind='linear')
                tckren[2] = scipy.interpolate.interp1d(ktab,rentab.T[1],kind='linear')
                tckren[4] = scipy.interpolate.interp1d(ktab,rentab.T[2],kind='linear')

                for nn in (0,2,4):
                    for a,b,Vab in itertools.izip(self.V[k][nn].row,self.V[k][nn].col,self.V[k][nn].data):
                        if a > b:
                            continue
                        elif a == b:
                            c = 1
                        else:
                            c = 2

                        Eab2= (self.basis[k][a].energy + self.basis[k][b].energy)/2.

                        coeff = tckren[nn](Eab2)
                        self.eigsrensubl[k][i] += c * coeff * cbar[a] * cbar[b] * Vab

        gc.collect()

    def vacuumE(self, ren="raw"):
        if ren=="raw":
            return self.eigenvalues[1][0]
        elif ren=="renlocal":
            return self.eigsrenlocal[1][0]
        elif ren=="rensubl":
            return self.eigsrensubl[1][0]
        else:
            raise ValueError("")
        # The vacuum is K-even

    def spectrum(self, k, ren="raw"):
        if ren=="raw":
            eigs = self.eigenvalues
        elif ren=="renlocal":
            eigs = self.eigsrenlocal
        elif ren=="rensubl":
            eigs = self.eigsrensubl
        else:
            raise ValueError("")

        if k==1:
            return scipy.array([x-self.vacuumE(ren=ren) for x in eigs[k][1:]])
        elif k==-1:
            return scipy.array([x-self.vacuumE(ren=ren) for x in eigs[k]])
        else:
            raise ValueError("")
            # Subtract vacuum energies
