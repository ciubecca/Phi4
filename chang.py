from scipy.optimize import fsolve
import numpy as np
from numpy import log, exp, sqrt, pi

def Xsol(x):
    """ Find X = g/M corresponding to x= g/m """
    func = lambda X1: 2/x**2+6/(pi*x)+1/X1**2-6/(pi*X1)+12*log(x/X1)/(pi**2)
    X0 = 6
    return fsolve(func, X0)[0]

def xsolmax(X):
    """ Return x dual to X in second branch """
    func = lambda x1: 2/x1**2+6/(pi*x1)+1/X**2-6/(pi*X)+12*log(x1/X)/(pi**2)
    x0 = 10
    return fsolve(func, x0)[0]

def xmintomax(x):
    """ Find corresponding x in second branch """
    return xsolmax(Xsol(x))

