#####################################################
#####################################################
#########                                   #########
#########   Integrands of VVV are defined   #########
#########                                   #########
#####################################################
#####################################################

import math
import scipy
import numpy
import vegas
from numpy import prod
from scipy.special import factorial
from VVVintegrands import *

pi = scipy.pi
m = 1.


#Heaviside theta
def HT(x):
    return 0.5 * (numpy.sign(x) + 1)

#Energy particle
def om(x):
    return math.sqrt(m**2+x**2)

#Symmetry factor
def sym(q,p,k):
    return (factorial(4)**3)/(factorial(q)**2)/(factorial(k)**2)/(factorial(p)**2)/factorial(4-q-k)/factorial(4-q-p)/factorial(4-p-k)




########################
######### Phi0 #########
########################

###################
##  diagram 0.1  ##
###################
def phi0_1(ET,eps,x):
    # I do a change of variables
    y= numpy.tan(x)
    jacobian = 1/(prod(numpy.cos(x))**2)
    #relativistic factor
    relfact= 1/((4*pi)**6)*1/(om(y[0])*om(y[1])*om(y[2])*om(y[3])*om(-y[0]-y[2]-y[3])*om(-y[1]-y[2]-y[3]))
    cut1= om(y[0])+om(-y[0]-y[2]-y[3])+om(y[2])+om(y[3])
    cut2= om(y[1])+om(-y[1]-y[2]-y[3])+om(y[2])+om(y[3])

    return sym(2,2,2)*relfact*jacobian*HT(cut1-ET)/(eps-cut1)*HT(cut2-ET)/(eps-cut2)

########################
######### Phi2 #########
########################

###################
##  diagram 2.1  ##
###################
def phi2_1(ET,eps,x):
    # I do a change of variables
    y= numpy.tan(x)
    jacobian = 1/(prod(numpy.cos(x))**2)
    #relativistic factor
    relfact= 1/((4*pi)**5)*1/(om(y[0])*om(y[1])*om(y[2])*om(-y[0]-y[2])*om(-y[1]-y[2]))
    cut1= om(y[0])+om(-y[0]-y[2])+om(y[2])
    cut2= om(y[1])+om(-y[1]-y[2])+om(y[2])

    return sym(2,2,1)*relfact*jacobian*HT(cut1-ET)/(eps-cut1)*HT(cut2-ET)/(eps-cut2)

###################
##  diagram 2.2  ##
###################
def phi2_2(ET,eps,x):
    # I do a change of variables
    y= numpy.tan(x)
    jacobian = 1/(prod(numpy.cos(x))**2)
    # (y[0],y[1],y[2])=(q1,k1,k2)
    relfact= 1/((4*pi)**5)*1/(om(y[0])*om(y[1])*om(y[2])*om(-y[1]-y[2])*om(-y[0]-y[1]-y[2]))

    cut1= om(y[1]+y[2])+om(y[1])+om(y[2])
    cut2= om(y[0])+om(y[0]+y[1]+y[2])+om(y[1])+om(y[2])

    return 2*sym(2,1,2)*relfact*jacobian*HT(cut1-ET)/(eps-cut1)*HT(cut2-ET)/(eps-cut2)

###################
##  diagram 2.3  ##
###################
def phi2_3(ET,eps,x):
    # I do a change of variables
    y= numpy.tan(x)
    jacobian = 1/(prod(numpy.cos(x))**2)
     # (y[0],y[1],y[2])=(k1,k2,k3)
    relfact= 1/((4*pi)**5)*1/(om(y[0])*om(y[1])*om(y[2])*om(-y[0]-y[1]-y[2])*om(-y[0]-y[1]-y[2]))

    cut1= om(y[0])+om(y[1])+om(y[2])+om(-y[0]-y[1]-y[2])
    cut2= cut1

    return sym(1,1,3)*relfact*jacobian*HT(cut1-ET)/(eps-cut1)*HT(cut2-ET)/(eps-cut2)

###################
##  diagram 2.4  ##
###################
def phi2_4(ET,eps,x):
    # I do a change of variables to map (-infinity,infinity) -> (-pi/2,pi/2)
    y= numpy.tan(x)
    jacobian = 1/(prod(numpy.cos(x))**2)
    # (y[0],y[1],y[2])=(p1,p2,p3)
    relfact= 1/((4*pi)**5)*1/(om(y[0])*om(y[1])*om(y[2])*om(-y[0]-y[1]-y[2])*om(y[0]+y[1]+y[2]))
    cut1= om(y[0])+om(y[1])+om(y[2])+om(-y[0]-y[1]-y[2])
    cut2= 2*om(-y[0]-y[1]-y[2])

    return 2*sym(1,3,1)*relfact*jacobian*HT(cut1-ET)/(eps-cut1)*HT(cut2-ET)/(eps-cut2)


########################
######### Phi4 #########
########################

###################
##  diagram 4.1  ##
###################
def phi4_1(ET,eps,x):
    # I do a change of variables
    y= numpy.tan(x)
    jacobian = 1/(prod(numpy.cos(x))**2)
    #relativistic factor
    relfact= 1/((4*pi)**4)*1/(om(y[0])*om(y[1])*om(-y[0])*om(-y[1]))
    cut1= om(y[0])+om(-y[0])
    cut2= om(y[1])+om(-y[1])

    return sym(2,2,0)*relfact*jacobian*HT(cut1-ET)/(eps-cut1)*HT(cut2-ET)/(eps-cut2)

###################
##  diagram 4.2  ##
###################
def phi4_2(ET,eps,x):
    # I do a change of variables
    y= numpy.tan(x)
    jacobian = 1/(prod(numpy.cos(x))**2)
    # (y[0],y[1])=(k1,q1)
    relfact= 1/((4*pi)**4)*1/(om(y[0])*om(y[1])*om(-y[0])*om(-y[1]))
    cut1= om(y[0])+om(-y[0])+om(-y[1])+om(y[1])
    cut2= om(y[0])+om(-y[0])

    return 2*sym(2,0,2)*relfact*jacobian*HT(cut1-ET)/(eps-cut1)*HT(cut2-ET)/(eps-cut2)

###################
##  diagram 4.3  ##
###################
def phi4_3(ET,eps,x):
    # I do a change of variables
    y= numpy.tan(x)
    jacobian = 1/(prod(numpy.cos(x))**2)
    # (y[0],y[1])=(k1,k2)
    relfact= 1/((4*pi)**4)*1/(om(y[0])*om(y[1])*om(y[0]+y[1])*om(y[0]+y[1]))
    cut1= om(y[0])+om(y[1])+om(-y[0]-y[1])
    cut2= cut1

    return sym(1,1,2)*relfact*jacobian*HT(cut1-ET)/(eps-cut1)*HT(cut2-ET)/(eps-cut2)

###################
##  diagram 4.4  ##
###################
def phi4_4(ET,eps,x):
    # I do a change of variables
    y= numpy.tan(x)
    jacobian = 1/(prod(numpy.cos(x))**2)
    # (y[0],y[1])=(k,p1)
    relfact= 1/((4*pi)**4)*1/(om(y[0])*om(y[1])*om(y[0])*om(-y[0]-y[1]))
    cut1= 2*om(y[0])
    cut2= om(y[0])+om(y[1])+om(y[0]+y[1])

    return 2*sym(2,1,1)*relfact*jacobian*HT(cut1-ET)/(eps-cut1)*HT(cut2-ET)/(eps-cut2)


###################
##  diagram 4.6  ##
###################
def phi4_6(ET,eps,x):
    # I do a change of variables
    y= numpy.tan(x)
    jacobian = 1/(prod(numpy.cos(x))**2)
    # (y[0],y[1])=(k1,k2)
    relfact= 1/((4*pi)**4)*1/(om(y[0])*om(y[1])*om(-y[0]-y[1])*m)
    cut1= om(y[0])+om(y[1])+om(-y[0]-y[1])
    cut2= cut1+m

    return 2*sym(0,1,3)*relfact*jacobian*HT(cut1-ET)/(eps-cut1)*HT(cut2-ET)/(eps-cut2)


########################
######### Phi6 #########
########################

###################
##  diagram 6.1  ##
###################
def phi6_1(ET,eps,x):
    # I do a change of variables
    y= numpy.tan(x)
    jacobian = 1/(prod(numpy.cos(x))**2)
    #y[0]=k2
    relfact= 1/((4*pi)**3)*1/(om(y[0])*om(y[0])*om(y[0]))
    cut1= 2*om(y[0])
    cut2= cut1

    return sym(1,1,1)*relfact*jacobian*HT(cut1-ET)/(eps-cut1)*HT(cut2-ET)/(eps-cut2)

###################
##  diagram 6.2  ##
###################
def phi6_2(ET,eps,x):
    # I do a change of variables
    y= numpy.tan(x)
    jacobian = 1/(prod(numpy.cos(x))**2)
    # y[0]=k1
    relfact= 1/((4*pi)**3)*1/(om(y[0])*om(-y[0])*m)
    cut1= 2*om(y[0])
    cut2= cut1+m

    return 2*sym(0,1,2)*relfact*jacobian*HT(cut1-ET)/(eps-cut1)*HT(cut2-ET)/(eps-cut2)

###################
##  diagram 4.5  ##
###################
# NOTE I changed this to a non-local correction
def phi0phi4_1(ET,eps,x):
    # I do a change of variables
    y= numpy.tan(x)
    jacobian = 1/(prod(numpy.cos(x))**2)
    # (y[0],y[1],y[2])=(k1,k2,k3)
    relfact= 1/((4*pi)**4)*1/(om(y[0])*om(y[1])*om(y[2])*om(-y[0]-y[1]-y[2]))
    cut1= om(y[0])+om(y[1])+om(y[2])+om(y[0]+y[1]+y[2])
    cut2= cut1

    return sym(0,0,4)*relfact*jacobian*HT(cut1-ET)/(eps-cut1)*HT(cut2-ET)/(eps-cut2)

############################
######### Phi2Phi4 #########
############################

#####################
##  diagram 2-4.1  ##
#####################
def phi2phi4_1(ET,eps,x):
    # I do a change of variables
    y= numpy.tan(x)
    jacobian = 1/(prod(numpy.cos(x))**2)
    #(y[0],y[1])=(k1,k2)
    relfact= 1/((4*pi)**3)*1/(om(y[0])*om(y[1])*om(-y[0]-y[1]))
    cut1= om(y[0])+om(y[1])+om(y[0]+y[1])
    cut2= cut1

    return sym(0,0,3)*relfact*jacobian*HT(cut1-ET)/(eps-cut1)*HT(cut2-ET)/(eps-cut2)

############################
######### Phi4Phi4 #########
############################

#####################
##  diagram 4-4.1  ##
#####################
def phi4phi4_1(ET,eps,x):
    # I do a change of variables
    y= numpy.tan(x)
    jacobian = 1/(prod(numpy.cos(x))**2)
    #(y[0],y[1])=(k1,k2)
    relfact= 1/((4*pi)**2)*1/(om(y[0])*om(y[-0]))
    cut1= 2*om(y[0])
    cut2= cut1

    return sym(0,0,2)*relfact*jacobian*HT(cut1-ET)/(eps-cut1)*HT(cut2-ET)/(eps-cut2)

