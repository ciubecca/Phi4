import sys
from scipy import integrate, prod, log, pi, array
import scipy
import vegas
from VVVintegrands import *
from scipy.special import factorial

#Heaviside theta
def HT(x):
    return 0.5 * (numpy.sign(x) + 1)

ET = 10
m = 1
eps = -1

neval = 10000
# XXX reset this to 10000
neval = 1000

def computeIntegral(nvar, integrand):
    cut=pi/2
    # integral variables mapped to the tangent to x:(-infinity,infinity)
    # replaced by x=tan(y) with y:(-pi/2,pi/2)
    nitn = 10

    ret = []
    integ = vegas.Integrator([[-cut,cut]]*nvar)

    # step 1 -- adapt to phi0_1; discard results
    integ(lambda x: integrand(ET,eps,x), nitn=nitn, neval=neval)
    result = integ(lambda x: integrand(ET,eps,x), nitn=nitn, neval=neval)

    return result.mean


print(computeIntegral(1, phi6_1))
print(computeIntegral(1, phi6_2))
