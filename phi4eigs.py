import inspect
import os
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
os.sys.path.insert(0,parentdir)

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

version = "v3_5-LV"
cutoff = 5.

def main(argv):
    if len(argv) < 3:
        print argv[0], "<fname> <Emax>"
        return -1
    
    fname = argv[1]
    Emax = float(argv[2])

    # Hardcoded parameters
    m=1.
    gmin = 0.
    gmax = 3.
    sigma = -30.
    npoints = 16
    neigs = num = 5

    a = phi1234.Phi1234()

    gList = scipy.linspace(gmin, gmax, num=npoints, endpoint=True)
    print gList

    a.loadMatrix(fname)

    a.buildBasis(k=1, Emax=Emax)
    a.buildBasis(k=-1, Emax=Emax)
    
    print 'K=1 full basis size = ',  a.fullBasis[1].size
    print 'K=-1 full basis size = ',  a.fullBasis[-1].size
    print 'K=1 basis size = ',  a.basis[1].size
    print 'K=-1 basis size = ',  a.basis[-1].size
    
    db = database.Database()

    for i, g in enumerate(gList):
        g = float(g)
        print "g = ", g

        existing = 0
        for k in (1,-1):
            for e in db.table.find(Emax=Emax, m=m, L=a.L, ren='rensubl', cutoff=cutoff, \
                                    k=k, dual=False, finiteVolCouplings=True, version=version):
                if (abs(e['g']-g)<10.**(-13.)) and (e['neigs']>=neigs):
                    existing += 1
                    break

        if existing >= 2:
            print 'Spectra for g=', g, ' already present'
            continue

        b = finiteVolH.FiniteVolH(a.L, m)
        g0, g2, g4 = b.directCouplings(g)
         
        a.setCouplings(g0=g0, g2=g2, g4=g4)
        print "Computing raw eigenvalues for g0,g2,g4 = ", a.g0,a.g2,a.g4

        a.computeHamiltonian(k=1, ren=False)
        a.computeHamiltonian(k=-1, ren=False)

        a.computeEigval(k=1, sigma=sigma, n=num, ren=False)
        a.computeEigval(k=-1, sigma=sigma, n=num, ren=False)
        
        print "Raw vacuum: ", a.vacuumE(ren="raw")
        
        a.renlocal(Er=a.vacuumE(ren="raw"))
        
        print "Computing renormalized eigenvalues for g0r,g2r,g4r = ", a.g0r,a.g2r,a.g4r
        
        a.computeHamiltonian(k=1, ren=True)
        a.computeHamiltonian(k=-1, ren=True)

        a.computeEigval(k=1, sigma=sigma, n=num, ren=True, corr=True, cutoff=cutoff)
        a.computeEigval(k=-1, sigma=sigma, n=num, ren=True, corr=True, cutoff=cutoff)

        print "Renlocal vacuum: ", a.vacuumE(ren="renlocal")
        print "Rensubl vacuum: ", a.vacuumE(ren="rensubl")

        for k in (1,-1):
            db.insert(k=k, Emax=Emax, L=a.L, ren="raw", g=g, spec=a.eigenvalues[k], basisSize=a.basis[k].size, \
                        neigs=neigs, version=version, dual=False, finiteVolCouplings=True, m=m)
            db.insert(k=k, Emax=Emax, L=a.L, ren="renlocal", g=g, spec=a.eigsrenlocal[k], basisSize=a.basis[k].size, \
                        neigs=neigs, Er=a.Er, cutoff=a.cutoff, version=version, dual=False, finiteVolCouplings=True, m=m)
            db.insert(k=k, Emax=Emax, L=a.L, ren="rensubl", g=g, spec=a.eigsrensubl[k], basisSize=a.basis[k].size, \
                        neigs=neigs, Er=a.Er, cutoff=a.cutoff, version=version, dual=False, finiteVolCouplings=True, m=m)

if __name__ == "__main__":
    main(sys.argv)
