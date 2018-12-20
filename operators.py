from profile_support import *
from operator import mul
from functools import reduce
import scipy
from math import factorial, floor, sqrt
import statefuncs
from statefuncs import Basis, Helper, tol
from collections import Counter
import itertools
from itertools import combinations, islice, permutations
from itertools import groupby
from scipy import exp, pi
from scipy.special import binom
import bisect
import numpy as np

debug = False

# XXX Review this function, to take Lambda into account?
def filterDlist(dlist, nd, ntot, helper):
    if nd==ntot:
        return tuple(sum([np.array(d) for d in dlist])) == (0,0)
    elif nd==ntot-1:
        return tuple(sum([np.array(d) for d in dlist])) in helper.allowedWn
    else:
        return True


def gendlists(state, nd, ntot, helper):
    """ Generates a list of all the possible combinations of momenta in the state that
    can be annihilated
    state: input state in representation 1
    nd: number of annihilation operators (number of modes to annihilate)
    ntot: total number of annihilation and creation operators
    allowedWn: all the allowed wave numbers in the basis
    """

    x = itertools.chain.from_iterable(([tuple(n)]*Zn for n,Zn in state))
    dlists = set(tuple(y) for y in combinations(x,nd))

    if debug and nd==0 and state ==[]:
        print("state 2", state)
        print("dlists 2", dlists)

    return (dlist for dlist in dlists if filterDlist(dlist, nd, ntot, helper))
