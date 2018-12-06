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

    def __init__(self, m, L, Emax, Lambda=np.inf):
        self.L = L
        self.m = m


        # Maximum sqrt(nx^2 + ny^2)
        nmax = floor(L/(2*pi)*min(Lambda, sqrt((Emax/2)**2-m**2)))
        self.nmax = nmax

        # self.omegaMat = array([[self._omega(array([nx,ny])) for nx in range(-nmax,nmax+1)] for ny in range(-nmax,nmax+1)])

        # Do not shift indices, but use negative indices for slight optimization, which avoids some operations in omega()
        self.omegaMat = np.zeros(shape=(2*self.nmax+1,2*self.nmax+1))
        for nx in range(-nmax,nmax+1):
            for ny in range(-nmax,nmax+1):
                self.omegaMat[nx][ny] = self._omega(array([nx,ny]))

        self.allowedWn = set()
        for nx in range(-nmax, nmax+1):
            for ny in range(-nmax, nmax+1):
                if sqrt(nx**2+ny**2)<=nmax:
                    self.allowedWn.add((nx,ny))


    # XXX This is slow
    def energy(self, state):
        """ Computes energy of state in Repr1 """
        return sum(Zn*self.omega(n) for n,Zn in state)

    def totwn(self, state):
        if state==[]:
            return array([0,0])
        return sum(Zn*n for n,Zn in state)

    def _omega(self, n):
        """ Energy corresponding to wavenumber n"""
        m = self.m
        return sqrt(m**2 + self.kSq(n))

    def omega(self, n):
        """ Energy corresponding to wavenumber n"""
        return self.omegaMat[n[0]][n[1]]
        # return self.omegaMat[n[0]+self.nmax][n[1]+self.nmax]

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

        self.helper = Helper(m, L, Emax, Lambda)
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
        allowedWn = helper.allowedWn

        ret = []

        for ny in itertools.count():
            for nx in itertools.count(1):

                n = array([nx,ny])

                if tuple(n) not in allowedWn:
                    if nx==1:
                        ret.sort(key=lambda n: omega(n))
                        return ret
                    else:
                        break

                ret.append(n)

        raise RuntimeError("Shouldn't get here")


    # @profile
    def _genNEstatelist(self, NEstate=[], idx=0):
        """ Recursive function generating all North-East moving states in Repr 1 starting from NEstate, by adding
        any number of particles with momentum self.NEwnlist[idx] """

        helper = self.helper
        m = helper.m
        Emax = self.Emax
        omega = helper.omega
        allowedWn = helper.allowedWn

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
                if tuple(WN) not in allowedWn or E+omega(WN)>Emax:
                    break
                newstate.append(mode)

            ret += self._genNEstatelist(newstate, idx+1)

        return ret

    def rotate(self, s):
        """ Rotate state counterclockwise by pi/2 """
        return [(np.dot(rot,n),Zn) for n,Zn in s]

    # @profile
    def buildBasis(self):
        """ Generates the basis starting from the list of RM states, in repr1 """

        helper = self.helper
        omega = helper.omega

        m = helper.m
        energy = helper.energy
        totwn = helper.totwn
        Emax = self.Emax
        allowedWn = helper.allowedWn
        occn = helper.occn


        # XXX Possible optimization: sort first by total wave number, and then by energy?

        NEsl1 = self.NEstatelist
        NEsl2 = [self.rotate(s) for s in NEsl1]
        NEsl3 = [self.rotate(s) for s in NEsl2]
        NEsl4 = [self.rotate(s) for s in NEsl3]

        # List of energies of the states in quadrants
        NEelist = [energy(s) for s in NEsl1]

        # Lists of total wavenumbers for states in quadrants
        NEwntotlist1 = [totwn(s) for s in NEsl1]
        NEwntotlist2 = [np.dot(rot,wntot) for wntot in NEwntotlist1]
        NEwntotlist3 = [np.dot(rot,wntot) for wntot in NEwntotlist2]
        NEwntotlist4 = [np.dot(rot,wntot) for wntot in NEwntotlist3]

        # Generate dictionary (kx,ky) -> ([idx]) , where idx is an index for the states in the South-East moving quadrant
        SEwnidx = dict()
        for idx, wn in enumerate(NEwntotlist4):
            if tuple(wn) not in SEwnidx.keys():
                SEwnidx[tuple(wn)] = [idx]
            else:
                SEwnidx[tuple(wn)].append(idx)

        ret = []

        # Rotate and join NE moving states counterclockwise
        for i1,s1 in enumerate(NEsl1):

            # Keep track of total energy
            E1 = NEelist[i1]
            # Keep track of total wave number
            WN1 = NEwntotlist1[i1]

            for i2,s2 in enumerate(NEsl2):
                E2 = E1 + NEelist[i2]
                WN2 = WN1 + NEwntotlist2[i2]

                # NEstatelist is ordered in energy
                if E2 > Emax:
                    break
                # We need to add at least another particle to have 0 total momentum.
                if tuple(WN2) not in allowedWn or E2+omega(WN2)>Emax:
                    continue

                # XXX Entering this inner cycle is the most expensive part
                for i3,s3 in enumerate(NEsl3):
                    E3 = E2 + NEelist[i3]
                    WN3 = WN2 + NEwntotlist3[i3]

                    # NEstatelist is ordered in energy
                    if E3>Emax:
                        break

                    # We need to add at least another particle to have 0 total momentum.
                    # Also, we cannot add anymore negative x momentum or positive y momentum in step 4
                    # XXX omega here is called many times, and it takes a long time in total
                    if WN3[0]>0 or WN3[1]<0 or tuple(WN3) not in allowedWn or E3+omega(WN3)>Emax:
                        continue

                    # There is no states that can cancel the total momentum
                    # XXX is this redundant?
                    if tuple(-WN3) not in SEwnidx.keys():
                        continue

                    for i4 in SEwnidx[tuple(-WN3)]:
                        E4 = E3 + NEelist[i4]
                        WN4 = WN3 + NEwntotlist4[i4]
                        s4 = NEsl4[i4]

                        if E4 > Emax:
                            break

                        if WN4[0]!=0 or WN4[1]!=0:
                            raise RuntimeError("Total momentum should be zero")

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
