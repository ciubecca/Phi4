import scipy
import scipy.sparse.linalg
import scipy.sparse
import math
import time
from operator import attrgetter
import gc
import numpy as np
import sys



##############################
##  new code for Phi4 basis and matrix generations
##  S.Rychkov Sep 2016
##############################

pi = scipy.pi


##############
# this code uses two state representations
# Representation 1:
# state is a list of the form:
#[(n,Zn),...] where n is a wavenumber and Zn occupation number
# NOTE: wavenumbers are always ordered in n
# NOTE: (n,Zn) is included only if Zn>0
#
# The basis is initially generated in Representation 1
##
# Representation 2:
# On the other hand for some subsequent manipulations it will be convenient to also use
# the representation of basis elements
# in the form of the throughout list:
# State =[ Zk, k=-nmax,...nmax ]
#
## Variables denoting state representations of the first( resp. second) type are always written starting
# from small (resp. capital) letter
#
# nmax is set later in the code after the basis is generated
#################

####
# utilitary functions for manipulating the states and computing their various characteristics
#####

def toState(state,nmax):
    # function transforming from the first to the second representation
    State = [0 for n in range(-nmax,nmax+1)]
    for (n,Zn) in state:
        State[nmax+n]=Zn
    return State

def omega(n,L,m):
    """ computes one particle energy from wavenumber"""
    return math.sqrt(m**2.+((2.*pi/L)*n)**2.)
    # check if this function consumes a lot of time, perhaps precompute it?

def k(n,L):
    """ computes momentum from wavenumber"""
    return (2.*pi/L)*n

def energy(state,L,m):
    en = 0
    for (n,Zn) in state:
        en+=Zn*omega(n,L,m)
    return en

def wavenumber(state):
    wn = 0
    for (n,Zn) in state:
        wn+=Zn*n
    return wn

def momentum(state,L):
    return (2*pi/L)*wavenumber(state)

def totocc(state):
    """ computes total occupation number """
    return sum([x[1] for x in state])

def reverseP(state):
    """ computes the image of the state under the momentum inversion """
    return [(-x[0],x[1]) for x in reversed(state)]

########################
##  combinatorial factor functions. These are slightly different from the old code
## (they take the list of occupation numbers and not the list of oscillator momenta)
## but achieve the same purpose
########################


def f4(occs):
    """ computes combinatorial factor for an ordered list of occupation numbers, total occ=4"""
    occs=sorted(occs)
    if occs==[1,1,1,1]:
        return 24
    elif occs==[1,1,2]:
        return 12
    elif occs==[2,2]:
        return 6
    elif occs==[1,3]:
        return 4
    else:
        return 1

def f3(occs):
    """ computes combinatorial factor for an ordered list of occupation numbers, total occ=3"""
    occs=sorted(occs)
    if occs==[1,1,1]:
        return 6
    elif occs==[1,2]:
        return 3
    else:
        return 1

def f2(occs):
    """ computes combinatorial factor for an ordered list of occupation numbers, total occ=2"""
    if occs==[1,1]:
        return 2
    else:
        return 1

#####################
### array containing products of factors for the Hamiltonian matrix
### coming from initial and final state (non-)invariance under P
##########

efactor = np.array([[1, math.sqrt(2.) ],[1/math.sqrt(2.) ,1]])
# call as
# efactor[pi,pj]
### i is that state on which one acts, j is the state one obtains
### here p = 1 if the state is P-inv, p=0 if not
### equivalent, slower, code is as follows :
# def efactors(pinvi,pinvj):
#     if pinvi==1:
#         extrafactor1=0.5
#     else:
#         extrafactor1=1./math.sqrt(2.)
#     if pinvj==1:
#         extrafactor2=2.
#     else:
#         extrafactor2=math.sqrt(2.)
#     return extrafactor1 * extrafactor2


class newPhi4():
    """ main FVHT class """
    def __init__(self):
        self.h0 = dict()
        # Free hamiltonian
        self.potential = dict()
        # Quartic operator (not multiplied by the coupling)

        self.eigenvalues = dict()
        # Dictionary of eigenvalues (for different K quantum numbers)
        self.eigenvectors = dict()

        self.totBasisSize = dict()
        # Basis size for each value of K

        self.compBasisSize = dict() #holds basis sizes used to actually compute the eigenvalues

        scipy.set_printoptions(precision=15)

        self.Z2eigs = (1,-1)

        # The two different K eigenvalues (for now we have only P=1)

#    @profile
    def buildBasisAndMatrix(self,L,ET,m):

        t0=time.clock()
        self.L=float(L)
        self.m=float(m)
        self.Emax=float(ET)
        self.ET=float(ET) #synonims


        ############################
        # Step 1: BUILDING BASIS
        # basis = list states at rest (momentum=0) with energy up to ET
        ############################
        # algorithm used is the same as in the old paper
        ############################################
        # Step 1a:  generate RMlist - list of right-moving (RM) states with energy up to ET, momentum up to ET/2
        # RMstates have strictly positive wavenumbers only
        # the RM vacuum is also included
        ##############################

        maxN1 = int(min((ET/2.)/k(1,L), ET/omega(1,L,m))) #maximal occupation number of n=1 mode

        RMlist = [[]] + [[(1,N)] for N in range(1,maxN1+1)]

        for n in range(2,10000): #go over all other modes, 10000 = infinity effectively
                imax = len(RMlist)
                stateadded = 0 #flag which controls termination of the cycle over nmax
                # if no states have been added then the flag remains zero and the cycle is terminated

                for i in range(imax): # cycle over all RMstates which already exist
                    RMstate = RMlist[i]
                    p0 = momentum(RMstate,L)
                    e0 = energy(RMstate,L,m)
                    maxNn = int(
                            math.floor(
                            min((ET/2. - p0)/k(n,L), (ET-e0)/omega(n,L,m))
                            )
                            )#maximal occupation number of mode n given the occupation numbers of all previous modes

                    for N in range(1,maxNn+1):
                        newstate=RMstate[:]
                        newstate.append((n,N)) #add all possible occupation numbers for mode n
                        RMlist.append(newstate)
                        stateadded = 1 # raise the flag that at least one state has been added

                if stateadded == 0:
                    break #no state has been added

        ########################
        # Step 1b: divide RMlist into a list of lists, RMdivided,
        # RMdivided[wn] contains all RM states of total wavenumber wn
        # also each RMdivided[wn] is ordered in energy
        #####################################

        nRMmax=max([wavenumber(RMstate) for RMstate in RMlist])
        RMdivided = [[] for ntot in range(nRMmax+1)] #initialize list of lists
        for RMstate in RMlist: #go over RMstates and append them to corresponding sublists
                RMdivided[wavenumber(RMstate)].append(RMstate)

        #now sort each sublist in energy
        for RMsublist in RMdivided:
                RMsublist.sort(key=lambda x: energy(x,L,m))


        ###########################
        # Step 1c: finally build the basis combining the RM, LM and zero mode states
        ##########################
        basis = {K:[] for K in (1,-1)} # K is the field parity eigenvalue
        ## i.e. K=1/-1 corresponds states with even/odd total number of particles

        self.basis = basis

        for nRM,RMsublist in enumerate(RMdivided):
                for i, RMstate in enumerate(RMsublist):
                    ERM = energy(RMstate,L,m)
                    occRM = totocc(RMstate)
                    for LMstate in RMsublist[i:]: # LM part of the state will come from the same sublist.
                        # we will just have to reverse it
                        # NOTA BENE: as before, we only generate states up to spatial reflection, and so
                        # take the position of LMState to be greater or equal to the position of RMstate

                        ELM = energy(LMstate,L,m)
                        deltaE = ET - ERM - ELM
                        if deltaE < 0: #if this happens, we can break since subsequent LMstates have even higher
                            #energy (RMsublist is ordered in energy)
                            break

                        maxN0 = int(math.floor(deltaE/m))
                        LMstate = reverseP(LMstate) #reverse momenta
                        occLM=totocc(LMstate)

                        for N0 in range(maxN0+1):
                            #possible values for the occupation value at rest
                            if N0>0:
                                state = LMstate + [(0,N0)] + RMstate
                            else:
                                state = LMstate + RMstate
                            occ = occLM + occRM + N0
                            if occ%2 == 0:
                                basis[1].append(state)
                            else:
                                basis[-1].append(state)

        ### sort basis in energy
        for K in (1,-1):
            basis[K].sort(key= lambda x: energy(x,L,m))

        #### compute maximal wavenumber
        nmax = max ([max([x[0] for x in state]) for state in basis[1][1:]+basis[-1]]) #exclude the vacuum

        ### compute dictionary of oscillator energies up to nmax, to reduce the number of arithmetic operators below
        Omega = {n: omega(n,L,m) for n in range(-nmax, nmax+1)}

        print "nmax=",nmax, ", Basis sizes: ",len(basis[1]),len(basis[-1])

        ### create some lists of state attributes to avoid repeated operations below
        energylist = {}
        Pinv = {}
        Basis = {}
        StatePos = {}
        for K in (1,-1):
            # which states are P-inv and which are not:
            Pinv[K] = [int(reverseP(state)==state) for state in basis[K]]

            # state energies:
            energylist[K] = [energy(state,L,m) for state in basis[K]]

            # list of states in Representation 2 (see the beginning of the file):
            Basis[K] = [toState(state,nmax) for state in basis[K]]

            # finally lookup dictionary:
            # notice a small improvement compared to the old work.
            # states and their P-reflections can be looked up in a single dictionary,
            # which returns the position of the state in basis
            preparedict = ([(toState(state,nmax),i) for i,state in enumerate(basis[K])] +
                           [(toState(reverseP(state),nmax),i) for i,state in enumerate(basis[K]) if Pinv[K][i]==0])

            StatePos[K] = {tuple(x[0]) : x[1] for x in preparedict}
            # NOTA BENE: only tuples are hashable
            # NOTA BENE: only states in Representation 2 will be looked up
            # notice that Representation 2 is generally longer (contains zeros) and thus hashing takes a bit longer
            # than if we were hashing Representation 1
            # but I checked that it's not a big overhead
            # on the other had Representation 2 will be crucial to have a manageable code for acting with oscillators


        ############################
        # Step 2: BUILDING PHI^4 operators
        ############################
        # we will need 4-creation operator, 3 creation+ 1 annihilation,
        # and (a part of) 2 creation + 2 annihilation parts of Phi^4
        # the missing parts can be obtained by transpositions

        ###################
        # Step 2a: 4-creation operator part
        # it's represented as a 4 particle state in Representation 1 format
        # we can just take the part of basis[1] which contains 4-particle states,
        # Since those states have been generated up to P-invariance, we have to add reversed states where needed
        ########################


        # Data structures of the matrix in the sparse COO format
        # Constructing the matrices in this form is cheap both memory-wise and cpu-wise
        data = {K:[] for K in (-1,1)}
        row = {K:[] for K in (-1,1)}
        col = {K:[] for K in (-1,1)}


        V40 = []
        for i,state in enumerate(basis[1]):
            if totocc(state)<>4:
                continue
            V40.append(state)

            if Pinv[1][i]==0:
                V40.append(reverseP(state))

        V40.sort(key= lambda x: energy(x,L,m)) # sort V40 in energy
        # precompute various lists to use below
        V40factor = [f4([x[1] for x in state])/(4*math.sqrt(np.prod([Omega[n]**Zn for (n,Zn) in state]))) for state in V40]
        V40energy = [energy(state,L,m) for state in V40]


        ######################
        # Step 2b: V31 part consisting of a^+(k1) a^+(k2) a^+(k3) a(k4), k1+k2+k3=k4.
        # We build it as a list of lists V31[k4] =[ state ]
        # where momentum(state)=k4, and state is in Representation 1
        # we only need to keep those operators which have energy(state)<=ET since others
        ##########################

        V31={k4: [] for k4 in range(-nmax,nmax+1)}
        for k4 in range(-nmax,nmax+1):
            # i,j,kk are k1,k2,k3
            for i in range(-nmax,nmax+1):
                for j in range(max(i,k4-i-nmax),
                               min(nmax,int(math.floor((k4-i)/2.)))+1):
                    kk = k4-i-j
                    # NOTE the boundaries for range of j ensure that j<=kk<=nmax

                    Emin = Omega[i]+Omega[j]+Omega[kk] #  The minimal energy of the state obtained by acting
                    if Emin > ET:
                        continue

                    # NB i<=j<=k
                    if i<>j and j<>kk:
                        V31[k4].append([(i,1),(j,1),(kk,1)])
                        continue
                    if i==j and j<>kk:
                        V31[k4].append([(i,2),(kk,1)])
                        continue
                    if i<>j and j==kk:
                        V31[k4].append([(i,1),(kk,2)])
                        continue
                    if i==kk:
                        V31[k4].append([(i,3)])

            V31[k4].sort(key= lambda x: energy(x,L,m))

        # precompute various lists to use below
        V31factor = {k4: [f3([x[1] for x in state])/(4*math.sqrt(Omega[k4]*np.prod([Omega[n]**Zn for (n,Zn) in state])))
                               for state in V31[k4]] for k4 in V31.keys()}
        V31energy = {k4: [energy(state,L,m)- Omega[k4] for state in V31[k4]] for k4 in V31.keys()}

        ########################
        # Step 2c: build the V22 part consisting of a^+(k1) a^+(k2) a(k3) a(k4), k1+k2=k3+k4.
        # We store it as a list of lists V22[k3,k4] =[ state ]
        # where momentum(state)=k3+k4 and state is in Representation 1 format
        # we only need to keep operators with energy(state)<= ET
        # we will only store lexicographically ordered part of V22, in the sense of (A.4) of our first paper
        # http://arxiv.org/abs/1412.3460
        # but also including diagonal part which will be separated below
        ########################

        V22={(k3,k4): [] for k3 in range(-nmax,nmax+1) for k4 in range(k3,nmax+1)}

        for k3 in range(-nmax,nmax+1):
            for k4 in range(k3,nmax+1):
                # i,j are k1,k2
                for i in range(max(-nmax,k3+k4-nmax),
                               min(nmax,int(math.floor((k3+k4)/2.)))+1):
                    j=k3+k4-i #range for i is such that i<=j<=nmax
                    Emin = Omega[i]+Omega[j] #the energy of the state will be at least this
                    if Emin> ET:
                        continue

                    if sorted([abs(i),abs(j)])>sorted([abs(k3),abs(k4)]):
                        continue #only consider lexicographically ordered part of V22, in the sense of (A.4)
                        # but also including diagonal part which will be separated below

                    if i<>j:
                        V22[(k3,k4)].append([(i,1),(j,1)])
                    else:
                        V22[(k3,k4)].append([(i,2)])

                V22[(k3,k4)].sort(key= lambda x: energy(x,L,m))

        def symfactor((k3,k4)):
            if k3==k4:
                 return 1
            else:
                 return 2

        V22factor = {key: [symfactor(key)*f2([x[1] for x in state])/
                                (4*math.sqrt(Omega[key[0]]*Omega[key[1]]*np.prod([Omega[n]**Zn for (n,Zn) in state])))
                                for state in V22[key]]
                                for key in V22.keys()}
        V22energy = {key: [energy(state,L,m) - Omega[key[0]] - Omega[key[1]] for state in V22[key]] for key in V22.keys()}

        ###############################
        # Step 3: BUILDING PHI^4 matrix in the basis
        ###############################

        # this is factor produced when action on state |occ0> by a^+(Zn)
        # saves time to isolate it
        # call as oscFactor[occ0][Zn]
        oscFactor = [[ math.sqrt(np.prod([occ+1 for occ in range(occ0,occ0+Zn)])) for Zn in range(5)]
                     for occ0 in range(0,int(ET))]


        ############################
        # Step 3a: act with V40 on every state
        ################################


        for K in (1,-1):
            for i,State in enumerate(Basis[K]):
                # it's very useful to have the initial state in Representation 2, since acting with any
                # oscilator is just an arithmetic operation like in line ***** below
                # doing the same in Representation 1 would be a mess
                en = energylist[K][i]

                for numv, v in enumerate(V40):
                    if en + V40energy[numv] <= ET:
                        newState = State[:]
                        factor = 1.

                        for (n,Zn) in v:
                            #occ0=newState[nmax+n]
                            factor *= oscFactor[newState[nmax+n]][Zn]
                            newState[nmax+n]+=Zn # *****

                        j = StatePos[K][tuple(newState)] # lookup the state

                        row[K].append(i)
                        col[K].append(j)
                        data[K].append(1./L*factor * V40factor[numv] * efactor[Pinv[K][i], Pinv[K][j]])

                        # NB important to accumulate (+=) since the same matrix element can be touched twice
                        # if j is not P-inv
                    else:
                        break #reached beyond the maximal energy

        ##########################
        # Step 3b: act with V31 on every state
        #########################

        for K in (1,-1):
            for i,State in enumerate(Basis[K]):
                en = energylist[K][i]
                for k4,Zk4 in basis[K][i]:
                    # NB: note how useful it is to have Representation 1 as well
                    # Indeed we can cycle now immediately only over the wavenumbers k4 for which Zk4 is nonzero
                    # and act with V31[k4] corresponding to these momenta
                    # if we had to cycle over all indices from -nmax to nmax to find which Zn are nonzero
                    # it would be a big overhead
                    for numv, v in enumerate(V31[k4]): # only acting with these elements will give nonzero results
                        if en + V31energy[k4][numv] <= ET:
                            newState = State[:]

                            factor = math.sqrt(newState[nmax+k4]) #factor produced by a(k4)

                            newState[nmax + k4] -= 1

                            for (n,Zn) in v:
                                #occ0=newState[nmax+n]
                                factor *= oscFactor[newState[nmax+n]][Zn]
                                newState[nmax+n] += Zn

                            j = StatePos[K][tuple(newState)]

                            row[K].append(i)
                            col[K].append(j)
                            data[K].append((4./L*factor * V31factor[k4][numv] * efactor[Pinv[K][i],Pinv[K][j]]))

                        else:
                            break #reached beyond the maximal energy for given initial state i


        ##########################
        # Step 3c: act with V22 on every state
        # we generate simultaneously diagonal and offdiagonal part
        #########################


        for K in (1,-1):
            for i,State in enumerate(Basis[K]):
                en = energylist[K][i]
                ks =[kk for kk,Zkk in basis[K][i]]
                for ii,k3 in enumerate(ks):
                    for k4 in ks[ii:]:
                        if k3==k4 and State[nmax+k3]<2:
                            continue

                        for numv, v in enumerate(V22[(k3,k4)]): # can only act with these elements
                            if en + V22energy[(k3,k4)][numv]  <= ET:
                                newState = State[:]
                                factor = math.sqrt(newState[nmax+k4]) #factor produced by a(k4)
                                newState[nmax+k4] -= 1
                                factor *= math.sqrt(newState[nmax+k3]) #factor produced by a(k3)
                                newState[nmax+k3] -= 1

                                for (n,Zn) in v:
                                    #occ0=newState[nmax+n]
                                    factor *= oscFactor[newState[nmax+n]][Zn]
                                    newState[nmax+n]+=Zn

                                j = StatePos[K][tuple(newState)]

                                row[K].append(i)
                                col[K].append(j)
                                data[K].append(6./L* factor * V22factor[(k3,k4)][numv]*efactor[Pinv[K][i],Pinv[K][j]])

                            else:
                                break #reached beyond the maximal energy

        self.potential = {}
        for K in (-1,1):
            # Construct the matrix in the COO format from the computed vectors
            M = scipy.sparse.coo_matrix((data[K],(row[K],col[K])),
                    shape=(len(basis[K]),len(basis[K])))
            # Add transpose and subtract diagonal part once.
            # In this step the elements of the matrix with the same row and column indices should be resummed
            M = M + M.transpose() - scipy.sparse.spdiags(M.diagonal(),0,len(basis[K]),len(basis[K]))


            self.potential[K] = scipy.sparse.csr_matrix(M)

        print "number of nonzero entries in the matrices:", self.potential[1].nnz, self.potential[-1].nnz
        print "time taken to construct matrices:", time.clock()-t0

        for K in (1,-1):
                self.totBasisSize[K] = len(self.basis[K])

        self.h0 = {K: energylist[K] for K in (1,-1)}


    def saveMatrix(self, fname):
            """ Saves the potential and free hamiltonian to file """
            t = (fname, self.L, self.Emax, self.m) + tuple([e for k in self.Z2eigs for e in (self.h0[k],self.potential[k].data,self.potential[k].row,self.potential[k].col)])
            scipy.savez(*t)

    def loadMatrix(self, fname):
            """ Loads the potential and free hamiltonian from file """
            f = scipy.load(fname)
            self.L = f['arr_0']
            self.Emax = f['arr_1']
            self.m = f['arr_2']

            for i, k in enumerate(self.Z2eigs):
                self.h0[k] = f['arr_'+str(3+i*4)]
                self.totBasisSize[k] = len(self.h0[k])
                self.potential[k] = scipy.sparse.coo_matrix((f['arr_'+(str(4+i*4))], (f['arr_'+(str(5+i*4))], f['arr_'+(str(6+i*4))])), shape=(self.totBasisSize[k], self.totBasisSize[k]))


    def computeEigval(self, k, g4, e=None, sigma=0, n=6):
            """ Sets the internal variables self.eigenvalues
            k: K-parity quantum number
            e: cutoff energy <= Emax. If None: e = Emax
            g4: quartic coupling, appearing in the Hamiltonian as  g4*\phi^4 (notice the normalization)
            n: number of eigenvalues to compute
            sigma: point around which we should look for the eigenvalues."""

            t0 = time.clock()
            if e == None:
                e = self.Emax

            if e > self.Emax:
                raise ValueError("Cutoff exceeds Emax")

            h0 = self.h0[k]
            potential = self.potential[k]
            basisSize = sum(1 for x in h0 if x <= e)
            # Dimension of the reduced basis with energy <= e

            hamiltonian = scipy.sparse.spdiags(h0[0:basisSize],0,basisSize,basisSize) + g4*potential.tocsr()[0:basisSize,0:basisSize]

            self.g4 = g4 #store for future reference
            self.compBasisSize[k]=basisSize

            v0 = scipy.zeros(basisSize)
            for i in range(10):
                v0[i]=1.

            (self.eigenvalues[k], eigenvectorstranspose) = scipy.sparse.linalg.eigsh(hamiltonian, n, v0=v0,
                                            which='SA', return_eigenvectors=True)

            self.eigenvectors[k] = eigenvectorstranspose.T
            # Relatively quick method to calculate the eigenvalues closest to sigma in absolute value, for a symmetric matrix.

            gc.collect()
            print "time elapsed:", time.clock()-t0

    def vacuumE(self):
            return self.eigenvalues[1][0]
            # The vacuum is K-even

    def spectrum(self, k):
            if k==1:
                return scipy.array([x-self.vacuumE() for x in self.eigenvalues[k][1:]])
            else:
                return scipy.array([x-self.vacuumE() for x in self.eigenvalues[k]])
            # Subtract vacuum energies

#a = newPhi4(); a.buildBasisAndMatrix(10. ,16. ,1.) #L,ET,m
