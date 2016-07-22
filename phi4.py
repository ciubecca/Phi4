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
from matrix import Matrix
from scipy import exp, pi


# I use this "tolerance" parameter throughout the code to
# avoid possible numerical issues when confronting energies with Emax
tol = 0.0001

def comb(*x):
    """ computes combinatorial factor for list of elements """
    return factorial(len(x))/scipy.prod(list(map(factorial,collections.Counter(x).values())))


class Phi4():
    """ main FVHT class """
    def __init__(self):
        self.L=None
        self.m=None

        self.basis = {}
        self.h0 = {}
        self.V = {k:{} for k in {-1,1}}

        self.eigenvalues = {"raw":{}, "renlocal":{}}
        self.eigenvectors = {"raw":{}, "renlocal":{}}

        self.compBasisSize = {1: None, -1:None}
        # Holds basis sizes used to actually compute the eigenvalues

        scipy.set_printoptions(precision=15)


    def buildBasis(self, k, L, m, Emax, occmax=None):
        """ Builds the full Hilbert space basis """
        self.L = float(L)
        self.m = float(m)

        self.basis[k] = Basis.fromScratch(L=self.L, Emax=Emax, m=self.m, k=k, occmax=occmax)

    def buildMatrix(self, k):
        """ Builds the full hamiltonian in the basis of the free hamiltonian.
        This is computationally intensive. It can be skipped by loading the matrix from file """
        L=self.L
        m=self.m

        basis = self.basis[k]
        lookupBasis = self.basis[k]
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

            self.V[k][n] = (offdiag_V+offdiag_V.transpose()+Matrix(lookupBasis, basis, diag_V)).to('coo')*self.L


    def setCouplings(self, g0, g2, g4):
        self.g0 = g0
        self.g2 = g2
        self.g4 = g4
        c = 2.*sum([1/(2.*pi)*scipy.special.kn(0,n*self.m*self.L) for n in range(1,10)])
        self.m1 = self.m*exp(-2.*pi*c)

    def renlocal(self, Emax, Er):
        self.g0r, self.g2r, self.g4r = \
            renorm.renlocal(self.g0, self.g2, self.g4, Emax, m=self.m1, Er=Er)

    def computeHamiltonian(self, k, Emax, ren, Er=0):
        if ren=="raw":
            V = self.V[k][0]*self.g0 + self.V[k][2]*self.g2 + self.V[k][4]*self.g4
        elif ren=="renlocal":
            V = self.V[k][0]*self.g0r + self.V[k][2]*self.g2r + self.V[k][4]*self.g4r
        else:
            raise ValueError()

        H0 = self.h0[k]
        H = H0 + V

        basisL = Basis.fromBasis(self.basis[k], lambda v: v.energy <= Emax)
        basisH = Basis.fromBasis(self.basis[k], lambda v: v.energy > Emax)

        Hll = H.sub(basisL, basisL)
        self.compBasisSize[k] = Hll.M.shape[0]

        # Choose "alpha" vectors
        basisAlpha = Basis.fromBasis(basisL, lambda v: any(v[0]==n and v.occ==n for n in (0,2,4)))

        # Construct new basis vectors
        propagator = (Er*Matrix(basisH, basisH, scipy.sparse.eye(basisH.size)) - H0.sub(basisH,basisH)).to("csc").inverse()
        psialpha1 = propagator*V.sub(basisH, basisAlpha)
        psialpha2 = V.sub(basisH, basisH)*psialpha1


        # Add tails
        # TODO add subleading tails
        # basisT = highBasis
        # self.Hlt =
        # self.Htt =

        return Hll.M

    def computeEigval(self, k, Emax, ren, sigma=0, neigs=10):
        """ Sets the internal variables self.eigenvalues
        k : field parity
        Emax : max energy of truncated Hilbert space
        ren : renormalization procedure "raw" or "renlocal"
        neigs : number of eigenvalues to compute
        sigma : point around which we should look for eigenvalues.
        """

        H = self.computeHamiltonian(k, Emax, ren)

        (self.eigenvalues[ren][k], eigenvectorstranspose) = scipy.sparse.linalg.eigsh(H, k=neigs, sigma=sigma,
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


    def saveMatrix(self, fname):
        """ Saves the potential and free hamiltonian to file """
        t = (fname, self.L, self.m, \
            self.fullBasis[1].Emax, self.fullBasis[1].nmax,  \
            self.fullBasis[-1].Emax, self.fullBasis[-1].nmax, \
            self.h0[1].M.data,self.h0[1].M.row,self.h0[1].M.col, \
            self.V[1][0].M.data,self.V[1][0].M.row,self.V[1][0].M.col, \
            self.V[1][2].M.data,self.V[1][2].M.row,self.V[1][2].M.col, \
            self.V[1][4].M.data,self.V[1][4].M.row,self.V[1][4].M.col, \
            self.h0[-1].M.data,self.h0[-1].M.row,self.h0[-1].M.col, \
            self.V[-1][0].M.data,self.V[-1][0].M.row,self.V[-1][0].M.col, \
            self.V[-1][2].M.data,self.V[-1][2].M.row,self.V[-1][2].M.col, \
            self.V[-1][4].M.data,self.V[-1][4].M.row,self.V[-1][4].M.col \
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

            self.h0[k] = Matrix(basisI, basisJ, scipy.sparse.coo_matrix((f['arr_'+(str(z+i*n))], (f['arr_'+(str(z+1+i*n))],
                f['arr_'+(str(z+2+i*n))])), shape=(basisI.size, basisJ.size)))
            self.V[k][0] = Matrix(basisI, basisJ, scipy.sparse.coo_matrix((f['arr_'+(str(z+3+i*n))],
                (f['arr_'+(str(z+4+i*n))], f['arr_'+(str(z+5+i*n))])), shape=(basisI.size, basisJ.size)))
            self.V[k][2] = Matrix(basisI, basisJ, scipy.sparse.coo_matrix((f['arr_'+(str(z+6+i*n))],
                (f['arr_'+(str(z+7+i*n))], f['arr_'+(str(z+8+i*n))])), shape=(basisI.size, basisJ.size)))
            self.V[k][4] = Matrix(basisI, basisJ, scipy.sparse.coo_matrix((f['arr_'+(str(z+9+i*n))],
                (f['arr_'+(str(z+10+i*n))], f['arr_'+(str(z+11+i*n))])), shape=(basisI.size, basisJ.size)))
