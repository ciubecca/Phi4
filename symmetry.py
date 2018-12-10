import numpy as np
from numpy import array, sqrt, e, pi

# Counter clockwise rotation by 90 degrees
rot = array([[0,-1],[1,0]])
# Counter clockwise rotation by 180 degrees
rot2 = array([[-1,0],[0,-1]])
# Counter clockwise rotation by 270 degrees
rot3 = array([[0,1],[-1,0]])
# Reflection wrt y axis
refly = array([[-1,0],[0,1]])
# Reflection wrt x axis
reflx = array([[1,0],[0,-1]])
# Reflection wrt diagonals
xs = array([0,1],[1,0])
ys = array([0,-1],[-1,0])

def rotate(s):
    """ Rotate state counterclockwise by pi/2 """
    return [(np.dot(rot,n),Zn) for n,Zn in s]

# Define all transformations
def I(s):
    return s
def S(s):
    return [(np.dot(rot,n),Zn) for n,Zn in s]
def S2(s):
    return [(np.dot(rot2,n),Zn) for n,Zn in s]
def S3(s):
    return [(np.dot(rot3,n),Zn) for n,Zn in s]
def X(s):
    return [(np.dot(refly,n),Zn) for n,Zn in s]
def Y(s):
    return [(np.dot(reflx,n),Zn) for n,Zn in s]
def XS(s):
    return [(np.dot(xs,n),Zn) for n,Zn in s]
def YS(s):
    return [(np.dot(ys,n),Zn) for n,Zn in s]

def toTuple(state):
    """ Transform state in Repr1 with numpy arrays to tuples to store momenta """
    return [(tuple(n), Zn) for n,Zn in state]

def genTransformed(state):
    """ Generate the set of all transformed states under all symmetries for any given Fock space state """
    return  {bytes(helper.torepr2(f(s))) for f in (I, S, S2, S3, X, Y, XS, YS)}



