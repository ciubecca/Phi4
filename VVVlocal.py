import math
import scipy
import numpy
import vegas
import time
import os
from numpy import prod
from scipy.special import factorial
from VVVintegrands import *


t0=time.clock()


pi = scipy.pi
m = 1.

length=10.

eps=-.1

#Heaviside theta
def HT(x):
    return 0.5 * (numpy.sign(x) + 1)

#Energy particle
def om(x):
    return math.sqrt(m**2+x**2)

#Symmetry factor
def sym(q,p,k):
    return (factorial(4)**3)/(factorial(q)**2)/(factorial(k)**2)/(factorial(p)**2)/factorial(4-q-k)/factorial(4-q-p)/factorial(4-p-k)


################################################
################################################
#########           Integrals          #########
################################################
################################################


cut=pi/2 #integrant variables maped to the tangent to x:(-infinity,finitny) replaced by x=tan(y) with y:(-pi/2,pi/2)

ETlist=numpy.arange(10, 20.5, 0.5)

nitn=10
neval=10000

ET=ETlist[0]


print("\n********************* Phi0 ************************** \n")

c0 = [0 for i in ETlist]
integ4 = vegas.Integrator([[-cut,cut], [-cut,cut], [-cut,cut], [-cut,cut]])

print("===> diagram 1 ...")
# step 1 -- adapt to phi0_1; discard results
integ4(lambda x: phi0_1(ET,eps,x), nitn=nitn, neval=neval)
for i,ET in enumerate(ETlist):
    # step 2 -- integ has adapted to phi0_1; keep results
    result = integ4(lambda x: phi0_1(ET,eps,x), nitn=nitn, neval=neval)
    #print(result.summary())
    #print('ET=%.1f result = %s chi2/dof=%.2f   Q = %.2f' % (ET, result, result.chi2/result.dof, result.Q))
    c0[i]+=result

for i,ET in enumerate(ETlist):
    print("ET=", ET," c0=", c0[i],)


print("\n********************* Phi2 ************************** \n")


c2 = [0 for i in ETlist]
integ3 = vegas.Integrator([[-cut,cut], [-cut,cut], [-cut,cut]])

print("===> diagram 1 ...")
# step 1 -- adapt to phiX_Y; discard results
integ3(lambda x: phi2_1(ET,eps,x), nitn=nitn, neval=neval)
for i,ET in enumerate(ETlist):
    # step 2 -- integ has adapted to phi0_1; keep results
    result = integ3(lambda x: phi2_1(ET,eps,x), nitn=nitn, neval=neval)
    #print('ET=%.1f result = %s chi2/dof=%.2f   Q = %.2f' % (ET, result, result.chi2/result.dof, result.Q))
    c2[i]+=result

## I reset the integrator settings before going to the next diagram
integ3.set()

print("===> diagram 2 ...")
integ3(lambda x: phi2_2(ET,eps,x), nitn=nitn, neval=neval)
for i,ET in enumerate(ETlist):
    result = integ3(lambda x: phi2_2(ET,eps,x), nitn=nitn, neval=neval)
    #print('ET=%.1f result = %s chi2/dof=%.2f   Q = %.2f' % (ET, result, result.chi2/result.dof, result.Q))
    c2[i]+=result

## I reset the integrator settings before going to the next diagram
integ3.set()

print("===> diagram 3 ...")
integ3(lambda x: phi2_3(ET,eps,x), nitn=nitn, neval=neval)
for i,ET in enumerate(ETlist):
    result = integ3(lambda x: phi2_3(ET,eps,x), nitn=nitn, neval=neval)
    #print('ET=%.1f result = %s chi2/dof=%.2f   Q = %.2f' % (ET, result, result.chi2/result.dof, result.Q))
    c2[i]+=result
## I reset the integrator settings before going to the next diagram
integ3.set()

print("===> diagram 4 ...")
integ3(lambda x: phi2_4(ET,eps,x), nitn=nitn, neval=neval)
for i,ET in enumerate(ETlist):
    result = integ3(lambda x: phi2_4(ET,eps,x), nitn=nitn, neval=neval)
    #print('ET=%.1f result = %s chi2/dof=%.2f   Q = %.2f' % (ET, result, result.chi2/result.dof, result.Q))
    c2[i]+=result

for i,ET in enumerate(ETlist):
    print("ET=", ET," c2=", c2[i],)



print("\n********************* Phi4 ************************** \n")


c4 = [0 for i in ETlist]

integ2 = vegas.Integrator([[-cut,cut], [-cut,cut]])

print("===> diagram 1 ...")
# step 1 -- adapt to phiX_Y; discard results
integ2(lambda x: phi4_1(ET,eps,x), nitn=nitn, neval=neval)
for i,ET in enumerate(ETlist):
    # step 2 -- integ has adapted to phi0_1; keep results
    result = integ2(lambda x: phi4_1(ET,eps,x), nitn=nitn, neval=neval)
    print('ET=%.1f result = %s chi2/dof=%.2f   Q = %.2f' % (ET, result, result.chi2/result.dof, result.Q))
    c4[i]+=result

## I reset the integrator settings before going to the next diagram
integ2.set()

print("===> diagram 2 ...")
integ2(lambda x: phi4_2(ET,eps,x), nitn=nitn, neval=neval)
for i,ET in enumerate(ETlist):
    result = integ2(lambda x: phi4_2(ET,eps,x), nitn=nitn, neval=neval)
    print('ET=%.1f result = %s chi2/dof=%.2f   Q = %.2f' % (ET, result, result.chi2/result.dof, result.Q))
    c4[i]+=result

## I reset the integrator settings before going to the next diagram
integ2.set()

print("===> diagram 3 ...")
integ2(lambda x: phi4_3(ET,eps,x), nitn=nitn, neval=neval)
for i,ET in enumerate(ETlist):
    result = integ2(lambda x: phi4_3(ET,eps,x), nitn=nitn, neval=neval)
    #print('ET=%.1f result = %s chi2/dof=%.2f   Q = %.2f' % (ET, result, result.chi2/result.dof, result.Q))
    c4[i]+=result
## I reset the integrator settings before going to the next diagram
integ2.set()

print("===> diagram 4 ...")
integ2(lambda x: phi4_4(ET,eps,x), nitn=nitn, neval=neval)
for i,ET in enumerate(ETlist):
    result = integ2(lambda x: phi4_4(ET,eps,x), nitn=nitn, neval=neval)
    #print('ET=%.1f result = %s chi2/dof=%.2f   Q = %.2f' % (ET, result, result.chi2/result.dof, result.Q))
    c4[i]+=result## I reset the integrator settings before going to the next diagram
integ2.set()

integ3tmp = vegas.Integrator([[-cut,cut], [-cut,cut], [-cut,cut]])
print("===> diagram 5 ...")
integ3tmp(lambda x: phi4_5(ET,eps,x), nitn=nitn, neval=neval)
for i,ET in enumerate(ETlist):
    result = integ3tmp(lambda x: phi4_5(ET,eps,x), nitn=nitn, neval=neval)
    #print('ET=%.1f result = %s chi2/dof=%.2f   Q = %.2f' % (ET, result, result.chi2/result.dof, result.Q))
    c4[i]+=result

print("===> diagram 6 ...")
integ2(lambda x: phi4_6(ET,eps,x), nitn=nitn, neval=neval)
for i,ET in enumerate(ETlist):
    result = integ2(lambda x: phi4_6(ET,eps,x), nitn=nitn, neval=neval)
    #print('ET=%.1f result = %s chi2/dof=%.2f   Q = %.2f' % (ET, result, result.chi2/result.dof, result.Q))
    c4[i]+=result

for i,ET in enumerate(ETlist):
    print("ET=", ET," c4=", c4[i],)



print("\n********************* Phi6 ************************** \n")

c6 = [0 for i in ETlist]
neval=1000

integ1 = vegas.Integrator([[-cut,cut]])

print("===> diagram 1 ...")
# step 1 -- adapt to phiX_Y; discard results
integ1(lambda x: phi6_1(ET,eps,x), nitn=nitn, neval=neval)
for i,ET in enumerate(ETlist):
    result = integ1(lambda x: phi6_1(ET,eps,x), nitn=nitn, neval=neval)
    print('ET=%.1f result = %s chi2/dof=%.2f   Q = %.2f' % (ET, result, result.chi2/result.dof, result.Q))
    c6[i]+=result

integ1.set()

print("===> diagram 2 ...")
integ1(lambda x: phi6_2(ET,eps,x), nitn=nitn, neval=neval)
for i,ET in enumerate(ETlist):
    result = integ1(lambda x: phi6_2(ET,eps,x), nitn=nitn, neval=neval)
    print('ET=%.1f result = %s chi2/dof=%.2f   Q = %.2f' % (ET, result, result.chi2/result.dof, result.Q))
    c6[i]+=result

for i,ET in enumerate(ETlist):
    print("ET=", ET," c6=", c6[i],)



print("\n******************* Phi2Phi4 ************************ \n")

c24 = [0 for i in ETlist]
neval=1000

integ24 = vegas.Integrator([[-cut,cut], [-cut,cut]])

print("===> diagram 1 ...")
integ24(lambda x: phi2phi4_1(ET,eps,x), nitn=nitn, neval=neval)
for i,ET in enumerate(ETlist):
    result = integ24(lambda x: phi2phi4_1(ET,eps,x), nitn=nitn, neval=neval)
    print('ET=%.1f result = %s chi2/dof=%.2f   Q = %.2f' % (ET, result, result.chi2/result.dof, result.Q))
    c24[i]+=result

for i,ET in enumerate(ETlist):
    print("ET=", ET," c24=", c24[i],)


print("\n******************* Phi4Phi4 ************************ \n")


neval=1000
c44 = [0 for i in ETlist]

integ44 = vegas.Integrator([[-cut,cut]])

print("===> diagram 1 ...")
integ44(lambda x: phi4phi4_1(ET,eps,x), nitn=nitn, neval=neval)
for i,ET in enumerate(ETlist):
    result = integ44(lambda x: phi4phi4_1(ET,eps,x), nitn=nitn, neval=neval)
    print('ET=%.1f result = %s chi2/dof=%.2f   Q = %.2f' % (ET, result, result.chi2/result.dof, result.Q))
    c44[i]+=result

for i,ET in enumerate(ETlist):
    print("ET=", ET," c44=", c44[i],)

t1=time.clock()

print(t1-t0)

