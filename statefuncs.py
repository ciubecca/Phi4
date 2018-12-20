import scipy
from scipy import array, pi, sqrt
from math import floor, factorial
import numpy as np
from symmetry import *
from helper import *


class Basis():
    """ Class used to store and compute a basis of states"""


    def __init__(self, k, stateList, helper, statePos, ncomp,
            repr1=True, repr1Emax=None):
        """ Standard constructor
        k: parity quantum number
        stateset: set or list of states in representation 1
        helper: Helper object
        ncomp: symmetry components of each state
        """
        self.k = k
        self.helper = helper
        self.Emax = helper.Emax
        self.Lambda = helper.Lambda
        self.repr1 = repr1
        self.size = len(stateList)

        if repr1==False:
            totwn = helper.totwn2
            energy = helper.energy2
            maxmom = helper.maxmom2
            occn = occn2
            self.repr1Emax = repr1Emax
        else:
            totwn = helper.totwn
            energy = helper.energy
            maxmom = helper.maxmom

        # Retrieve the transformation of the indices to get lists sorted in energy
        energyList = [energy(state) for state in stateList]
        idx = np.argsort(np.array(energyList))
        # Reverse indices
        idx2 = {j:i for i,j in enumerate(idx)}
        # Remap the indices
        self.stateList = [stateList[idx[i]] for i in range(self.size)]
        self.energyList = [energyList[idx[i]] for i in range(self.size)]
        self.statePos = {state: idx2[i] for state,i in statePos.items()}

        # Symmetry types
        self.ncomp = [ncomp[idx[i]] for i in range(self.size)]

        self.occnList = [occn(state) for state in self.stateList]
        # Maximal single particle momentum in each state
        self.maxmom = [maxmom(s) for s in self.stateList]

        # Check assumptions
        el = self.energyList
        assert  all(el[i] <= el[i+1]+tol for i in range(len(el)-1))
        assert (max(el) <= self.Emax+tol)
        assert (max(self.maxmom) <= self.Lambda+tol)
        assert all(sum(totwn(s)**2)==0 for s in self.stateList)
        assert all(1-2*(occn(state)%2)==k for state in self.stateList)

        # Convert some states to representation 1
        if self.repr1==False and self.repr1Emax!=None:
            self.stateList2 = self.stateList[:]
            self.stateList = []
            for i,e in enumerate(self.stateList2):
                if e>Emax+tol:
                    break
                self.stateList.append(helper.torepr1(self.stateList2[i]))



    def subidxlist(self, Emax, Lambda=np.inf):
        """ Return the indices of states within smaller cutoffs """
        return [i for i in range(self.size) if self.energyList[i]<=Emax+tol and self.maxmom[i]<=Lambda+tol]


    @classmethod
    def fromScratch(self, m, L, Emax, Lambda=np.inf):
        """ Builds the truncated Hilbert space up to cutoff Emax from scratch, in repr1
        m: mass
        L: side of the torus
        Emax: maximal energy of the states
        """

        self.helper = Helper(m, L, Emax, Lambda)
        helper = self.helper
        m = helper.m
        energy = helper.energy

        self._occmax = int(floor(Emax/m)+tol)

        self._buildBasis(self)

        # Make the representation of each state unique by sorting the oscillators
        self.bases = {k: [sortOsc(s) for s in self.bases[k]] for k in (-1,1)}

        return {k: self(k, self.bases[k], helper, self.statePos[k], ncomp=self.ncomp[k]) for k in (-1,1)}


    def __len__(self):
        return len(self.stateList)

    def __repr__(self):
        return str([toCanonical(s) for s in self.stateList])

    def _genNEwnlist(self, Emax, Lambda):
        """ Generate list of North-East moving wave numbers momenta, nx > ny >= 0,
        sorted in energy, below cutoffs Emax and Lambda """

        helper = self.helper
        omega = helper.omega
        allowedWn = helper.allowedWn

        self.NEwnlist = []

        for ny in itertools.count():
            for nx in itertools.count(1):

                n = array([nx,ny])

                if tuple(n) not in allowedWn:

                    if nx==1:
                        self.NEwnlist.sort(key=lambda n: omega(n))
                        return
                    else:
                        break

                self.NEwnlist.append(n)

        raise RuntimeError("Shouldn't get here")


    def _genNEstatelist(self, NEstate=[], idx=0):
        """ Recursive function generating all North-East moving states in Repr 1 starting from NEstate, by adding
        any number of particles with momentum self.NEwnlist[idx] """

        helper = self.helper
        m = helper.m
        Emax = helper.Emax
        omega = helper.omega
        allowedWn = helper.allowedWn
        energy = helper.energy
        totwn = helper.totwn
        minEnergy = helper.minEnergy

        if idx == len(self.NEwnlist):
            # TODO Could Sort oscillators inside state here
            return [NEstate]

        # Two-dimensional NE-moving wave number
        n = self.NEwnlist[idx]

        # Keep track of current energy
        E = energy(NEstate)
        # Keep track of current total wave number
        WN = totwn(NEstate)

        ret = []

        for Zn in itertools.count():
            newstate = NEstate[:]

            if Zn>0:
                mode = (n, Zn)
                E += omega(n)
                WN += n
                # We need to add at least another particle to have 0 total momentum.
                # XXX Check
                # if tuple(WN) not in allowedWn or E+minEnergy(WN)>Emax+tol:
                if E+minEnergy(WN)>Emax+tol:
                    break
                newstate.append(mode)

            ret += self._genNEstatelist(self, newstate, idx+1)

        return ret


    def _buildBasis(self):
        """ Generates the basis starting from the list of RM states, in repr1 """

        helper = self.helper
        omega = helper.omega

        m = helper.m
        energy = helper.energy
        totwn = helper.totwn
        Emax = helper.Emax
        Lambda = helper.Lambda
        allowedWn = helper.allowedWn
        minEnergy = helper.minEnergy

        # Generate list of all NE moving momenta
        self._genNEwnlist(self, Emax, Lambda)

        # Generate list of all NE moving states, and sort them by energy
        NEstatelist = self._genNEstatelist(self)
        NEstatelist.sort(key=lambda s: energy(s))
        # List of energies of the states in quadrants
        NEelist = [energy(s) for s in NEstatelist]

        NEsl1 = NEstatelist
        NEsl2 = [rotate(s) for s in NEsl1]

        # Lists of total wavenumbers for states in quadrants
        NEwntotlist1 = [totwn(s) for s in NEsl1]
        NEwntotlist2 = [np.dot(rot,wntot) for wntot in NEwntotlist1]

        # North-moving states in first and second quadrants
        Nsl12 = {}

        # Rotate and join NE moving states counterclockwise
        for i1,s1 in enumerate(NEsl1):

            # Keep track of total energy
            E1 = NEelist[i1]
            # Keep track of total wave number
            WN1 = NEwntotlist1[i1]

            # XXX The commented part is wrong. Can we save any more time?
            # for i2 in range(i1, len(NEsl2)):
                # s2 = NEsl2[i2]
            for i2,s2 in enumerate(NEsl2):

                E2 = E1 + NEelist[i2]
                WN2 = WN1 + NEwntotlist2[i2]

                # NEstatelist is ordered in energy
                if E2 > Emax+tol:
                    break

                # We need to add at least another particle to have 0 total momentum.
                # XXX CHECK
                if E2+minEnergy(WN2)>Emax+tol:
                    continue

                s12 = s1+s2

                if tuple(WN2) not in Nsl12.keys():
                    Nsl12[tuple(WN2)] = [s12]
                else:
                    Nsl12[tuple(WN2)].append(s12)

        # Create states moving north or south, and join them pairwise
        Nsl12 = [list(sorted(states, key=energy)) for states in Nsl12.values()]
        N12elist = [[energy(s) for s in states] for states in Nsl12]
        N12occlist = [[occn(s) for s in states] for states in Nsl12]

        # List of states in Representation 1, which are not related by symmetries
        self.bases = {k:[] for k in (-1,1)}
        idx = {k:0 for k in (-1,1)}
        # Dictionary of indices for states in Representation 2, including the "redundant" ones by symmetry
        self.statePos = {k:{} for k in (-1,1)}
        # Number of Fock states in the singlet representation of state. This is used to compute the appropriate normalization
        # factors
        self.ncomp = {k:[] for k in (-1,1)}

        for i, states in enumerate(Nsl12):

            for j1, s in enumerate(states):
                s34 = rotate(rotate(s))
                e34 = N12elist[i][j1]
                o34 = N12occlist[i][j1]

                # Save some time by using reflection
                for j2 in range(j1, len(states)):
                    s12 = states[j2]
                    e = e34 + N12elist[i][j2]
                    o = o34 + N12occlist[i][j2]

                    if e > Emax+tol:
                        break

                    for Z0 in range(int(floor((Emax-e)/m+tol))+1):
                        occtot = o + Z0
                        k = 1-2*(occtot%2)

                        if Z0==0:
                            state = s34 + s12
                            # The state already exists (for every Z0), when taking symmetries into account
                            if bytes(helper.torepr2(state)) in self.statePos[k]:
                                break
                        else:
                            state = s34 + s12 + [(array([0,0]),Z0)]

                        transStates = genTransformed(state, helper)

                        # Number of Fock space states in the singlet state
                        self.ncomp[k].append(len(transStates))
                        self.bases[k].append(toCanonical(state))

                        for s in transStates:
                            self.statePos[k][s] = idx[k]
                        idx[k] += 1

        return
