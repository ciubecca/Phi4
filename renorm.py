import sys
from scipy import integrate, prod, log, pi, array
import scipy
import vegas
from VVVintegrands import *
from scipy.special import factorial

def ft0(g2, g4, E, m):
    return (g2**2./pi + g4**2.*(-3/(2*pi) + 18/pi**3 * log(E/m)**2))/(E**2.)


def ft2(g2, g4, E, m):
    return (g2*g4*12./pi + g4**2.*72./pi**2. * log(E/m))/(E**2.)


def ft4(g2, g4, E, m):
    return (g4**2. * 36./pi) / (E**2)


#Heaviside theta
def HT(x):
    return 0.5 * (numpy.sign(x) + 1)


#Symmetry factor
def sym(q,p,k):
    return (factorial(4)**3)/(factorial(q)**2)/(factorial(k)**2)/(factorial(p)**2)/factorial(4-q-k)/factorial(4-q-p)/factorial(4-p-k)


class renVV2():

    def __init__(self, g4, EL, eps, m=1, g2=0):

        self.VV2 = {}
        self.EL = EL
        self.m = m
        self.eps = eps

        self.VV2[0] = integrate.quad(lambda E: ft0(g2,g4,E,m)/(eps-E),EL,scipy.inf)[0]
        self.VV2[2] = integrate.quad(lambda E: ft2(g2,g4,E,m)/(eps-E),EL,scipy.inf)[0]
        self.VV2[4] = integrate.quad(lambda E: ft4(g2,g4,E,m)/(eps-E),EL,scipy.inf)[0]


    def om(x):
        return sqrt(self.m**2+x**2)

class renVV3Test():
    def __init__(self, glist):
        self.VV3loc = {n : {g:1 for g in glist} for n in (0,2,4,6)}
        self.V0V4 = {g: 1 for g in glist}
        self.V2V4 = {g: 1 for g in glist}
        self.V4V4 = {g: 1 for g in glist}



class renVV3():
    def __init__(self, glist, m, eps, ET):
        self.neval = 10000

        self.glist = glist

        self.m = m
        self.eps = eps

        self.VV3loc = {}

# Local correction to identity
        self.VV3loc[0] = self.computeIntegral(4, phi0_1)

# Local correction to mass operator
        self.VV3loc[2] = self.computeIntegral(3, phi2_1)
        self.VV3loc[2] += self.computeIntegral(3, phi2_2)
        self.VV3loc[2] += self.computeIntegral(3, phi2_3)
        self.VV3loc[2] += self.computeIntegral(3, phi2_4)

# Local corrections to V4 NOTE I moved diagram 4.5 to the non-local corrections
        self.VV3loc[4] = self.computeIntegral(2, phi4_1)
        self.VV3loc[4] += self.computeIntegral(2, phi4_2)
        self.VV3loc[4] += self.computeIntegral(2, phi4_3)
        self.VV3loc[4] += self.computeIntegral(2, phi4_4)
        self.VV3loc[4] += self.computeIntegral(2, phi4_6)


        self.neval = 1000
# Local corrections to V6
        self.VV3loc[6] = self.computeIntegral(1, phi6_1)
# NOTE This is overestimated
        self.VV3loc[6] += self.computeIntegral(1, phi6_2)


# Convert to dict
        self.VV3loc = {n: {g: self.VV3loc[n][i] for i,g in enumerate(self.glist)}
            for n in self.VV3loc.keys()}

        self.neval = 10000
        # Bilocal corrections
        self.V0V4 = self.computeIntegral(3, phi0phi4_1)
        self.V2V4 = self.computeIntegral(2, phi2phi4_1)
        self.V4V4 = self.computeIntegral(1, phi4phi4_1)

        # Convert to dict
        self.V0V4 = {ET: self.V0V4[i] for i,g in enumerate(self.glist)}
        self.V2V4 = {ET: self.V2V4[i] for i,g in enumerate(self.glist)}
        self.V4V4 = {ET: self.V4V4[i] for i,g in enumerate(self.glist)}


    def computeIntegral(self, nvar, integrand):
        cut=pi/2
        # integral variables mapped to the tangent to x:(-infinity,infinity)
        # replaced by x=tan(y) with y:(-pi/2,pi/2)
        neval = self.neval
        nitn = 10
        eps = self.eps[glist[0]]
        ET = self.ET

        ret = []
        integ = vegas.Integrator([[-cut,cut]]*nvar)

        # step 1 -- adapt to phi0_1; discard results
        integ(lambda x: integrand(ET,eps,x), nitn=nitn, neval=neval)

        for g in self.glist:
        # step 2 -- integ has adapted to phi0_1; keep results
            result = integ(lambda x: integrand(ET,eps[g],x), nitn=nitn, neval=neval)
            ret.append(result.mean)

        return array(ret)


    def om(x):
        return sqrt(self.m**2+x**2)
