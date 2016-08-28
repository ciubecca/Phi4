import scipy
from scipy import pi, sqrt, array
from statefuncs import omega, State


pi=scipy.pi # this is pi=3.14159...

# NOT IMPORTANT. I use this "tolerance" parameter throughout the code to avoid possible numerical issues when confronting energies with Emax
tol = 0.0001

""" P denotes spatial parity, while K field parity. For now only the P-even sector is implemented """

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
            n*=state[i]
            state[i]-=1

        for i in self.clist:
            n*=state[i]+1
            state[i]+=1

        return (n, state)

    def apply(self, v, lookupbasis):
        """ Takes a state index in basis, returns another state index (if it
        belongs to the lookupbasis) and a proportionality coefficient.
        Otherwise raises LookupError
        lookupbasis can be different from basis,
        but it's assumed that they have the same nmax"""

        # if self.deltaE+v.energy < 0.-tol or self.deltaE+v.energy > lookupbasis.Emax+tol:
            # # FIX (ADD CHECK ON OCCMAX)
            # # The trasformed element surely does not belong to the basis if E>Emax
            # raise LookupError()

        n, newstate = self._transformState(v)
        # if n==0:
            # return (0, None)
        j = lookupbasis.lookup(newstate)

        c = 1.
        # Required for state normalization
        if v.isParityEigenstate(): c /= sqrt(2.)
        if lookupbasis[j].isParityEigenstate(): c *= sqrt(2)

        return (c*sqrt(n)*self.coeff, j)
