import scipy
import math
from scipy import integrate
from scipy import optimize
log = scipy.log
pi = scipy.pi
sqrt = scipy.sqrt
exp = scipy.exp

cutoff = 20.

class FiniteVolH:
    def __init__(self,L,m):
        self.L = L
        self.m = m

    def finiteVolCouplings(self, mx, g0, g2, g4):
        # mx can be either the original mass or the dual (infinite-volume) mass

        def f(x, a):
            return 1./(sqrt(a**2.+x**2.)*(exp(sqrt(a**2.+x**2.))-1.))

        L = self.L
        z = 1./pi*integrate.quad(lambda x: f(x,L*mx),0.,cutoff)[0]
        E0 = -1./(pi*L)*integrate.quad(lambda x: x**2.*f(x,L*mx),0.,cutoff)[0]
   
        return (g0 + E0/L + g2*z + 3.*g4*z**2.,
                g2 + 6*g4*z, # This does not include the m^2/2 term in H0
                g4)

 
    def directCouplings(self, g4):
        return self.finiteVolCouplings(self.m, 0, 0, g4)


    def dualCouplings(self, g4):

        def Z(mu):
            return 1/(4.*pi)*log(self.m**2./mu**2.)
        def Y(mu):
            return 1/(8.*pi)*(mu**2.-self.m**2.)   
        def f(mu):
            return self.m**2./2. + 6.*g4*Z(mu) + mu**2./4

        mu0 = g4/self.m
        self.mu = max(scipy.optimize.fsolve(f, mu0))

        Y = Y(self.mu)
        Z = Z(self.mu)

        g0 = 1/2.*self.m**2.*Z + 3.*g4*Z**2. + Y
        g2 = -0.75*self.mu**2.
        
        return self.finiteVolCouplings(self.mu, g0, g2, g4)
