import finiteVolH
import phi1234
import renorm
import sys
import scipy
import math
import database

def main(argv):
    if len(argv) < 2:
        print(argv[0], " <Emax> <g>")
        return -1

    Emax = float(argv[1])
    g = float(argv[2])

    # Hardcoded parameters
    m=1.
    sigma = -30.
    neigs = 1
    Llist = scipy.linspace(5, 16, 23)
    print(Llist)

    db = database.Database()

    for L in Llist:
        a = phi1234.Phi1234()

        a.buildFullBasis(k=1,Emax=Emax,m=m,L=L)
        # a.buildFullBasis(k=-1,Emax=Emax,m=m,L=L)
        print("k=1 basis size :", a.fullBasis[1].size)

        a.buildMatrix(k=1)
        a.buildBasis(k=1, Emax=Emax)

        approxQuery = {"g":g, "L":L, "Emax":Emax}
        exactQuery = {"k":1}
        if db.getObjList('spec', approxQuery=approxQuery, exactQuery=exactQuery) != []:
            print("Eigenvalues already present")
            continue

        b = finiteVolH.FiniteVolH(a.L, m)
        g0, g2, g4 = b.directCouplings(g)

        a.setCouplings(g0=g0, g2=g2, g4=g4)
        print("Computing raw eigenvalues for g0,g2,g4 = ", a.g0,a.g2,a.g4)

        a.computeHamiltonian(k=1, ren="raw")
        # a.computeHamiltonian(k=-1, ren=False)

        a.computeEigval(k=1, sigma=sigma, n=neigs, ren="raw")
        # a.computeEigval(k=-1, sigma=sigma, n=neigs, ren=False)

        print("Raw vacuum: ", a.vacuumE(ren="raw"))

        a.renlocal(Er=a.vacuumE(ren="raw"))

        print("Computing renormalized eigenvalues for g0r,g2r,g4r = ", a.g0r,a.g2r,a.g4r)

        a.computeHamiltonian(k=1, ren="renlocal")
        # a.computeHamiltonian(k=-1, ren=True)

        a.computeEigval(k=1, sigma=sigma, n=neigs, ren="renlocal")
        # a.computeEigval(k=-1, sigma=sigma, n=neigs, ren=True, corr=True, cutoff=cutoff)

        print("Renlocal vacuum: ", a.vacuumE(ren="renlocal"))
        # print "Rensubl vacuum: ", a.vacuumE(ren="rensubl")

        for k in (1,):
            for ren in ("raw","renlocal"):
                db.insert(k=k, Emax=Emax, L=a.L, ren=ren, g=g, spec=a.eigenvalues[ren][k],
                        eigv=a.eigenvectors[ren][k], basisSize=a.basis[k].size, neigs=neigs)

if __name__ == "__main__":
    main(sys.argv)
