import math
import scipy
import vegas
import numpy as np
from numpy import prod, sqrt
from scipy.special import factorial

pi = np.pi
m = 1

def sym(q,p,k):
    """ Symmetry factor """
    return (factorial(4)**3)/(factorial(q)**2)/(factorial(k)**2)/(factorial(p)**2)/factorial(4-q-k)/factorial(4-q-p)/factorial(4-p-k)


def HT(x):
    """ Heaviside theta """
    return 0.5 * (np.sign(x) + 1.)

def om(x,y):
    """ Energy particle """
    return sqrt(m**2+x**2+y**2)

def P2(x,y):
    """ momentum squared """
    return sqrt(x**2+y**2)


class Integrator():

    def __init__(self, nitn=20, neval=50000):
        self.nitn = nitn
        self.neval = neval

    def do(self, lam):
        cut = pi/2
        integ = vegas.Integrator([[-cut,cut], [-cut,cut], [-cut,cut], [-cut,cut], [-cut,cut], [-cut,cut]])

        # step 1 -- adapt to integrand; discard results
        integ(lambda x: self.integrand(x, lam), nitn=self.nitn, neval=self.neval)

        # step 2 -- integ has adapted to phi0_1; keep results
        ret = integ(lambda x: self.integrand(x, lam), nitn=self.nitn, neval=self.neval)
        return ret



class Phi0_1(Integrator):
    """ O(VV) vacuum diagram """
    sym = -1/24
    # Partial relativistic factor
    relfact = 1/pow(2*pi,6)*1/pow(2,4)

    def __init__(self, *args, **kwargs):
        super(Phi0_1, self).__init__(*args, **kwargs)


    def integrand(self, x, lam):
        """ O(VV) vacuum diagram:
            x: vector of momenta
            lam: momentum cutoff """

        sym = self.sym
        # Change of variables
        y = np.tan(x)
        jacobian = 1/(prod(np.cos(x))**2)
        # Relativistic factor
        relfact = self.relfact*1/(om(y[0],y[1])*om(y[2],y[3])*om(y[4],y[5])*om(-y[0]-y[2]-y[4],-y[1]-y[3]-y[5]))
        cut1 = om(y[0],y[1])+om(y[2],y[3])+om(y[4],y[5])+om(-y[0]-y[2]-y[4],-y[1]-y[3]-y[5])

        return sym * relfact * jacobian * 1/(cut1) *\
                HT(lam-P2(y[0],y[1])) *\
                HT(lam-P2(y[2],y[3])) *\
                HT(lam-P2(y[4],y[5])) *\
                HT(lam-P2(-y[0]-y[2]-y[4],-y[1]-y[3]-y[5]))

    def counterterm(self, lam):
        return 0
