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


def filterDlist(dlist, nd, ntot, helper):
    """
    Returns True if the given tuple of oscillators can be assigned to annihilation operators for the given type operator
    It takes into account, for instance, the maximal energy of the basis, or whether the momenta sum up to zero
    dlist: tuple of annihilation momenta
    nd: number of annihilation operators
    ntot: total number of operators in the composite operator
    helper: Helper object
    """
    if nd==ntot:
        return tuple(sum([np.array(d) for d in dlist])) == (0,0)
    elif nd==ntot-1:
        return tuple(sum([np.array(d) for d in dlist])) in helper.allowedWn
    else:
        return True


def gendlists(state, nd, ntot, helper):
    """ Generates a list of all the possible combinations of momenta in the state that
    can be annihilated
    nd: number of annihilation operators (number of modes to annihilate)
    state: input state in representation 1
    ntot: total number of annihilation and creation operators
    helper: Helper object
    """

    x = itertools.chain.from_iterable(([tuple(n)]*Zn for n,Zn in state))
    dlists = set(tuple(y) for y in combinations(x,nd))

    return (dlist for dlist in dlists if filterDlist(dlist, nd, ntot, helper))
