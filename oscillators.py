import scipy
from scipy import pi, sqrt, array
from statefuncs import omega, State


class NormalOrderedOperator():
    """ abstract class for normal ordered operator """
    def __init__(self,clist,dlist,L,m,nmax,extracoeff=1):
        self.clist = clist
        self.dlist = dlist
        self.L = L
        self.m = m
        self.coeff = extracoeff/\
            scipy.product([sqrt(2.*L*omega(n,L,m)) for n in clist+dlist])
        self.deltaE = sum([omega(n,L,m) for n in clist])-sum([omega(n,L,m) for n in dlist])
        # self.deltaN = len(clist)-len(dlist)

    def __repr__(self):
        return str(self.clist)+" "+str(self.dlist)


    def _transformState(self, state0):
        state = State(state0.occs[:], state0.nmax, fast=True)
        n=1.
        for i in self.dlist:
            if state[i]==0:
                raise LookupError()
                #return(0,None)
            n *= state[i]
            state[i] -= 1

        for i in self.clist:
            n *= state[i]+1
            state[i] += 1

        return (n, state)

    def apply(self, v, lookupbasis):
        """ Takes a state index in basis, returns another state index (if it
        belongs to the lookupbasis) and a proportionality coefficient.
        lookupbasis can be different from basis,
        but it's assumed that they have the same nmax"""

        n, newstate = self._transformState(v)
        j = lookupbasis.lookup(newstate)

        c = 1.
        # Required for state normalization
        if v.isParityEigenstate(): c /= sqrt(2.)
        if lookupbasis[j].isParityEigenstate(): c *= sqrt(2)

        return (c*sqrt(n)*self.coeff, j)
