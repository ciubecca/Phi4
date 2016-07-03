import finiteVolH
import phi1234
import renorm
import sys
import scipy
import math
from scipy import optimize
import json
import cProfile
import database

cutoff = 5.

def main(argv):
    if len(argv) < 3:
        print argv[0], " <Emax> <g> <fname>"
        return -1

    Emax = float(argv[1])
    g = float(argv[2])
    fname = argv[3]

    # Hardcoded parameters
    m=1.
    sigma = -30.
    neigs = 5

    a = phi1234.Phi1234()

    a.loadMatrix(fname)
    L = a.L

    a.buildBasis(k=1, Emax=Emax)
    a.buildBasis(k=-1, Emax=Emax)

    print 'K=1 basis size = ',  a.basis[1].size
    print 'K=-1 basis size = ',  a.basis[-1].size

    db = database.Database()

    for e in db.table.find(Emax=Emax, m=m, L=a.L, cutoff=cutoff):
        if (abs(e['g']-g)<10.**(-13.)) and (e['neigs']>=neigs):
            print 'Spectra for g=', g, ' already present'

    b = finiteVolH.FiniteVolH(a.L, m)
    g0, g2, g4 = b.directCouplings(g)

    a.setCouplings(g0=g0, g2=g2, g4=g4)
    print "Computing raw eigenvalues for g0,g2,g4 = ", a.g0,a.g2,a.g4

    a.computeHamiltonian(k=1, ren=False)
    a.computeHamiltonian(k=-1, ren=False)

    a.computeEigval(k=1, sigma=sigma, n=neigs, ren=False)
    a.computeEigval(k=-1, sigma=sigma, n=neigs, ren=False)

    print "Raw vacuum: ", a.vacuumE(ren="raw")

    # a.renlocal(Er=a.vacuumE(ren="raw"))

    # print "Computing renormalized eigenvalues for g0r,g2r,g4r = ", a.g0r,a.g2r,a.g4r

    # a.computeHamiltonian(k=1, ren=True)
    # a.computeHamiltonian(k=-1, ren=True)

    # a.computeEigval(k=1, sigma=sigma, n=num, ren=True, corr=True, cutoff=cutoff)
    # a.computeEigval(k=-1, sigma=sigma, n=num, ren=True, corr=True, cutoff=cutoff)

    # print "Renlocal vacuum: ", a.vacuumE(ren="renlocal")
    # print "Rensubl vacuum: ", a.vacuumE(ren="rensubl")

    for k in (1,-1):
        db.insert(k=k, Emax=Emax, L=a.L, ren="raw", g=g, spec=a.eigenvalues[k], basisSize=a.basis[k].size, neigs=neigs, m=m)
        # db.insert(k=k, Emax=Emax, L=a.L, ren="renlocal", g=g, spec=a.eigsrenlocal[k], basisSize=a.basis[k].size, \
                        # neigs=neigs, Er=a.Er, cutoff=a.cutoff, m=m)
        # db.insert(k=k, Emax=Emax, L=a.L, ren="rensubl", g=g, spec=a.eigsrensubl[k], basisSize=a.basis[k].size, \
                        # neigs=neigs, Er=a.Er, cutoff=a.cutoff, m=m)

if __name__ == "__main__":
    main(sys.argv)
