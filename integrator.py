import math
import scipy
import vegas
import numpy as np
from numpy import prod, sqrt, cos, sin
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

def omk(r):
    """ Energy for particle with radial momentum r """
    return sqrt(m**2 + r**2)

def P2(x,y):
    """ momentum squared """
    return sqrt(x**2+y**2)


class Integrator():

    def __init__(self, nitn=20, neval=50000):
        self.nitn = nitn
        self.neval = neval

    def do(self, lam):
        cut = pi/2
        integ = vegas.Integrator([[0,lam], [0,lam], [0,lam], [0,2*pi], [0,2*pi]])

        # step 1 -- adapt to integrand; discard results
        integ(lambda x: self.integrand(x, lam), nitn=self.nitn, neval=self.neval)

        # step 2 -- integ has adapted to phi0_1; keep results
        ret = integ(lambda x: self.integrand(x, lam), nitn=self.nitn, neval=self.neval)

        # Old version by Joan
        # integ = vegas.Integrator([[-cut,cut], [-cut,cut], [-cut,cut], [-cut,cut], [-cut,cut], [-cut,cut]])
        # integ(lambda x: self.integrandJoan(x, lam), nitn=self.nitn, neval=self.neval)
        # ret = integ(lambda x: self.integrandJoan(x, lam), nitn=self.nitn, neval=self.neval)

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
            lam: momentum cutoff
            the variable s is [arctan(r0), arctan(r1), arctan(r2), theta1, theta2]
            """

        sym = self.sym
        th = x[3:]

        # Change of variables
        # r = np.tan(x[:3])
        # jacobian = 1/(cos(x[0])*cos(x[1])*cos(x[2]))**2
        r = x[:3]
        jacobian = r[0]*r[1]*r[2]

        # Radial momentum of 4th particle
        r3 = sqrt((r[0]+r[1]*cos(th[0])+r[2]*cos(th[1]))**2 + (r[1]*sin(th[0])+r[2]*sin(th[1]))**2)
        # Energy of four particles
        e0 = omk(r[0])
        e1 = omk(r[1])
        e2 = omk(r[2])
        e3 = omk(r3)

        # Relativistic factor
        relfact = self.relfact * 1/(e0*e1*e2*e3)
        cut = e0+e1+e2+e3

        # 2 pi is to account for the omitted angle
        return (2*pi) * sym * relfact * jacobian * 1/(cut) * HT(lam-r3)

    def integrandJoan(self, x, lam):
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
