import sys
import matplotlib.pyplot as plt
import scipy
import math
import json
from scipy import pi, log, array, sqrt
from math import factorial
from matplotlib import rc
import database
from statefuncs import Basis
import itertools

occmax = 5
n = 3
output = "pdf"

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)


""" Takes the a coefficient of the discretized wafe functions and renormalizes it
    according to the basis element.
    Works for general particle number """
# XXX possible mistakes here
def normalizeWF(c, v):
    # Bose symmetry
    c *= scipy.prod([factorial(n) for n in v.occs])/factorial(v.occ)
    # Normalization of Fock space states
    c *= 1/sqrt(scipy.prod([factorial(n) for n in v.occs]))
    # Spatial parity symmetry
    if v.isParityEigenstate() == False:
        c *= 1/sqrt(2)

    return c

def main(argv):
    args = " <g> <Emax> <L>"
    if len(argv) < 4:
        print(argv[0], args)
        return -1

    g = float(argv[1])
    Emax = float(argv[2])
    L = float(argv[3])

    plt.figure(1)
    params = {'legend.fontsize': 8}
    plt.rcParams.update(params)

    db = database.Database()
    exactQuery = {"ren":"raw", "k":-1}
    approxQuery = {"g":g, "Emax":Emax, "L":L}
    eigv = db.getObjList("eigv", exactQuery=exactQuery, approxQuery=approxQuery)[0]

    basis = Basis.fromScratch(m=1, L=L, Emax=Emax, k=-1, occmax=occmax)

    wf = []
    basisnp = []
    for i,v in enumerate(basis):
        # Select only n particles basis states
        if v.occ==n:
            basisnp.append(v)
            wf.append(eigv[0][i])

    # Renormalize wave function
    wf = array([normalizeWF(c,v) for c,v in zip(wf, basisnp)])

    # Construct variables in the form [k1,...,kn,f(k1,...,kn)]
    data = []
    for c, v in zip(wf, basisnp):
        # List of momenta of wavenumbers in the state
        wavenumbers = itertools.chain(*[[wn]*v[wn] for wn in v.wnList()])

        # Take all possible inequivalent (n-1)-tuples of wave numbers, including symmetrization
        s = set(map(lambda x: itertools.islice(x,n-1), itertools.permutations(wavenumbers)))
        # Take also the opposite wavenumbers
        s |= set(map(,s))

        for a,b in s:
            data.append(array([a*2*pi/L,b*2*pi/L,c]))

    data = array(data)
    scipy.savetxt("data/wf{0:d}p_g={1:g}_L={2:g}_E={3:g}.csv".format(n,g,L,Emax),
            data.reshape(1,data.size),delimiter=",")


if __name__ == "__main__":
    main(sys.argv)
