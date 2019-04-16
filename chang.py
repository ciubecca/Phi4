# File containing various utility functions related to Chang duality


from scipy.optimize import fsolve
import numpy as np
from numpy import log, exp, sqrt, pi
from scipy.optimize import brentq

# Approximate self-dual point
gstar = 2.76194


def Xsol(x):
    """ Find X = g/M corresponding to x= g/m """
    func = lambda X1: 2/x**2+6/(pi*x)+1/X1**2-6/(pi*X1)+12*log(x/X1)/(pi**2)
    X0 = 6
    return fsolve(func, X0)[0]


def xsolmax(X):
    """ Return x dual to X in second branch """
    func = lambda x1: 2/x1**2+6/(pi*x1)+1/X**2-6/(pi*X)+12*log(x1/X)/(pi**2)
    xmax = 20
    # return max(fsolve(func, x0))
    return brentq(func, a=gstar, b=xmax)


def xmintomax(x):
    """ Find corresponding x in second branch """
    return xsolmax(Xsol(x))


def factorToSecondBranch(x):
    """ Proportional factor by which I have to multiply mass gap in first
    branch in order to get prediction in second branch via Chang duality """
    return xmintomax(x)/x


def xmintomax2(x):
    """ Direct solution, find effective smaller g corresponding to g in first branch """
    func = lambda y: 1/x**2+3/(pi*x)-1/y**2-3/(pi*y)+6*log(x/y)/(pi**2)
    xmax = 20
    return brentq(func, a=gstar, b=xmax)


def factorToSecondBranch2(x):
    """ Proportional factor by which I have to multiply mass gap in first
    branch in order to get prediction in second branch via Chang duality """
    return xmintomax2(x)/x

