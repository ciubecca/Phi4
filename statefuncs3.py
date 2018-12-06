import bisect
from sys import getsizeof as sizeof
import scipy
from scipy import array, pi, sqrt
from math import floor, factorial
from collections import Counter
from itertools import combinations
import itertools
import numpy as np

# Counter clockwise rotation by 90 degrees
rot = array([[0,-1],[1,0]])

class Helper():
    """ This is just a "helper" class used to conveniently compute energies of
    oscillators and states and so on"""

    def __init__(self, m, L):
        self.L = L
        self.m = m

    def energy(self, state):
        """ Computes energy of state in Repr1 """
        return sum(Zn*self.omega(n) for n,Zn in state)

    def totwn(self, state):
        if state==[]:
            return array([0,0])
        return sum([Zn*n for n,Zn in state])

    def omega(self, n):
        m = self.m
        return sqrt(m**2 + self.kSq(n))

    def kSq(self, n):
        """ Squared momentum corresponding to wave number n """
        L = self.L
        return (2*pi/L)**2*np.sum(n**2)

    def occn(self, s):
        """ Occupation number of state """
        return sum([Zn for n,Zn in s])


class Basis():
    """ Class used to store and compute a basis of states"""

    def __init__(self, m, L, k, Emax, Lambda=np.inf):
        """ Builds the truncated Hilbert space up to cutoff Emax from scratch, in repr1
        m: mass
        L: side of the torus
        Emax: maximal energy of the states
        """

        self.helper = Helper(m, L)
        helper = self.helper
        energy = helper.energy
        totwn = helper.totwn
        occn = helper.occn
        self.Emax = Emax
        m = helper.m
        self.k = k

        self.NEwnlist = self._genNEwnlist(Emax, Lambda)

        self.NEstatelist = self._genNEstatelist()
        self.NEstatelist.sort(key=lambda s: energy(s))

        self.stateList = self.buildBasis()
        self.elist = [energy(s) for s in self.stateList]

        # Check assumptions
        el = self.elist
        assert  all(el[i] <= el[i+1] for i in range(len(el)-1))
        assert (max(el) <= Emax)
        assert all(sum(totwn(s)**2)==0 for s in self.stateList)
        assert all(1-2*(occn(state)%2)==k for state in self.stateList)


    def __len__(self):
        return len(self.stateList)

    def reprState(self, state):
        return [(tuple(n), Zn) for n,Zn in state]

    def __repr__(self):
        return str([self.reprState(s) for s in self.stateList])

    def _genNEwnlist(self, Emax, Lambda):
        """ Generate list of North-East moving wave numbers momenta, nx > ny >= 0,
        sorted in energy, below cutoffs Emax and Lambda """

        helper = self.helper
        omega = helper.omega
        kSq = helper.kSq

        ret = []

        for ny in itertools.count():
            for nx in itertools.count(1):

                n = array([nx,ny])

                if omega(n)>Emax/2 or kSq(n)>Lambda**2:
                    if nx==1:
                        ret.sort(key=lambda n: omega(n))
                        return ret
                    else:
                        break

                ret.append(n)

        raise RuntimeError("Shouldn't get here")


    def _genNEstatelist(self, NEstate=[], idx=0):
        """ Recursive function generating all North-East moving states in Repr 1 starting from NEstate, by adding
        any number of particles with momentum self.NEwnlist[idx] """

        helper = self.helper
        m = helper.m
        Emax = self.Emax
        omega = helper.omega

        if idx == len(self.NEwnlist):
            # TODO Sort oscillators inside state !
            return [NEstate]

        # Two-dimensional NE-moving wave number
        n = self.NEwnlist[idx]

        # Keep track of current energy
        E = helper.energy(NEstate)
        # Keep track of current total wave number
        WN = helper.totwn(NEstate)

        ret = []

        for Zn in itertools.count():
            newstate = NEstate[:]

            if Zn>0:
                mode = (n, Zn)
                E += omega(n)
                WN += n
                # We need to add at least another particle to have 0 total momentum.
                if E + omega(WN) > Emax:
                    break
                newstate.append(mode)

            ret += self._genNEstatelist(newstate, idx+1)

        return ret

    def rotate(self, s):
        """ Rotate state counterclockwise by pi/2 """
        return [(np.dot(rot,n),Zn) for n,Zn in s]

    def buildBasis(self):
        """ Generates the basis starting from the list of RM states, in repr1 """

        helper = self.helper
        omega = helper.omega
        m = helper.m
        energy = helper.energy
        totwn = helper.totwn
        Emax = self.Emax
        occn = helper.occn

        NEsl1 = self.NEstatelist
        NEsl2 = [self.rotate(s) for s in NEsl1]
        NEsl3 = [self.rotate(s) for s in NEsl2]
        NEsl4 = [self.rotate(s) for s in NEsl3]

        ret = []

        # Rotate and join NE moving states counterclockwise
        for s1 in NEsl1:

            # Keep track of total energy
            E1 = energy(s1)
            # Keep track of total wave number
            WN1 = totwn(s1)

            for s2 in NEsl2:
                E2 = E1+energy(s2)
                WN2 = WN1+totwn(s2)

                # NEstatelist is ordered in energy
                if E2 > Emax:
                    break
                # We need to add at least another particle to have 0 total momentum.
                if E2+omega(WN2) > Emax:
                    continue

                for s3 in NEsl3:
                    E3 = E2+energy(s3)
                    WN3 = WN2+totwn(s3)

                    # NEstatelist is ordered in energy
                    if E3 > Emax:
                        break
                    # We need to add at least another particle to have 0 total momentum. Also, we cannot add anymore
                    # negative x momentum in step 4
                    if E3+omega(WN3)>Emax or WN3[0]>0:
                        continue

                    for s4 in NEsl4:
                        E4 = E3+energy(s4)
                        WN4 = WN3+totwn(s4)

                        if E4 > Emax:
                            break
                        # TODO We could use a dict in this step to select only states with appropriate momentum,
                        # instead of cycling through all of them
                        if (WN4**2).sum() != 0:
                            continue

                        # Add zero modes
                        for Z0 in itertools.count():
                            Etot = E4 + Z0*m
                            if Etot > Emax:
                                break

                            state = s1+s2+s3+s4
                            if Z0>0:
                                state += [(array([0,0]),Z0)]

                            if self.k == 1-2*(occn(state)%2):
                                ret.append(state)

        return sorted(ret, key=energy)
