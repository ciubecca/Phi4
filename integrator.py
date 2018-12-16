import math
import scipy
import vegas
import numpy as np
from numpy import prod, sqrt, cos, sin, log
from scipy.special import factorial

pi = np.pi
m = 1

def HT(x):
    """ Heaviside theta """
    return 0.5 * (np.sign(x) + 1.)

def om(r):
    """ Energy for particle with radial momentum r """
    return sqrt(m**2 + r**2)


class Integrator():

    def __init__(self, nitn=20, neval=50000):
        self.nitn = nitn
        self.neval = neval

    def do(self, lam):

        integ = vegas.Integrator(self.interval(lam))
        # step 1 -- adapt to integrand; discard results
        integ(lambda x: self.integrand(x, lam), nitn=self.nitn, neval=self.neval)
        # step 2 -- integ has adapted; keep results
        ret = integ(lambda x: self.integrand(x, lam), nitn=self.nitn/2, neval=self.neval*2)

        return ret



class Phi0_1(Integrator):
    """ O(g^2) vacuum diagram """
    # Total factor. g is normalized so it is divided by 4! in the Lagrangian
    factor = 1/factorial(4)*1/(2**4*(2*pi)**6)

    def __init__(self, *args, **kwargs):
        super(Phi0_1, self).__init__(*args, **kwargs)

    def interval(self, lam):
        return [[0,lam], [0,lam], [0,lam], [0,2*pi], [0,2*pi]]

    def integrand(self, x, lam):
        """ x: vector of momenta
            lam: momentum cutoff
            the variable s is [arctan(r0), arctan(r1), arctan(r2), theta1, theta2]
            """
        th = x[3:]
        # Change of variables
        r = x[:3]
        jacobian = r[0]*r[1]*r[2]
        # Radial momentum of 4th particle
        r3 = sqrt((r[0]+r[1]*cos(th[0])+r[2]*cos(th[1]))**2 + (r[1]*sin(th[0])+r[2]*sin(th[1]))**2)
        # Energy of four particles
        e0 = om(r[0])
        e1 = om(r[1])
        e2 = om(r[2])
        e3 = om(r3)
        # Non-relativistic propagator
        prop = -1/(e0+e1+e2+e3)
        # 2 pi is to account for the omitted angle
        return (2*pi) * self.factor * jacobian * prop * HT(lam-r3) * 1/(e0*e1*e2*e3)

    def counterterm(self, lam):
        # XXX This is wrong
        return -1/(48*(4*pi)**3)*lam



class Phi0_2(Integrator):
    """ O(g^2) correction to overlap """

    # There is an additional factor 1/2 compared to the energy correction
    factor = 1/factorial(4)*1/(2**4*(2*pi)**6)*1/2

    def __init__(self, *args, **kwargs):
        super(Phi0_2, self).__init__(*args, **kwargs)

    def interval(self, lam):
        return [[0,lam], [0,lam], [0,lam], [0,2*pi], [0,2*pi]]

    def integrand(self, x, lam):
        th = x[3:]
        r = x[:3]
        jacobian = r[0]*r[1]*r[2]
        r3 = sqrt((r[0]+r[1]*cos(th[0])+r[2]*cos(th[1]))**2 + (r[1]*sin(th[0])+r[2]*sin(th[1]))**2)
        e0 = om(r[0])
        e1 = om(r[1])
        e2 = om(r[2])
        e3 = om(r3)
        prop = 1/(e0+e1+e2+e3)
        # The propagator is squared
        return -(2*pi) * self.factor * jacobian * prop**2 * HT(lam-r3) * 1/(e0*e1*e2*e3)


class Phi1_1(Integrator):
    """ O(g^2) correction to g2.  The Hamiltonian is normalized as
    g2 V2 + g4/(4!) V4 """
    factor = 1/(4*factorial(4))*1/((2*pi)**4)

    def __init__(self, *args, **kwargs):
        super(Phi1_1, self).__init__(*args, **kwargs)

    def interval(self, lam):
        return [[0,lam], [0,lam], [0,2*pi]]

    def integrand(self, x, lam):
        th = x[2]
        r = x[:2]
        jacobian = r[0]*r[1]
        # Radial momentum of 3rd particle
        r2 = sqrt((r[0]+r[1]*cos(th))**2 + (r[1]*sin(th))**2)
        e0 = om(r[0])
        e1 = om(r[1])
        e2 = om(r2)
        # There are two diagrams in OFPT
        prop = 1/(m-(e0+e1+e2)) + 1/(-m-(e0+e1+e2))
        # 2 pi is to account for the omitted angle
        return (2*pi) * self.factor * jacobian * prop * HT(lam-r2) * 1/(e0*e1*e2)

    def counterterm(self, lam):
        # XXX This is wrong
        return -1/(12*(4*pi)**2)*log(lam/m)

