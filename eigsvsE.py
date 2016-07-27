import finiteVolH
import phi4
import renorm
import sys
import scipy
import math
import database

addTails = True
saveondb = True
json = False

# Hardcoded parameters
m = 1.
Emin = 5
Elim = 20
sigma = -30.
neigs = 1
k = 1

def main(argv):
    if len(argv) < 5:
        print(argv[0], "<L> <Emaxbar> <g> <k> <occmax>")
        return -1

    L = float(argv[1])
    Emaxbar = float(argv[2])
    g = float(argv[3])
    k = int(argv[4])
    try:
        occmax = int(argv[5])
    except IndexError:
        occmax = None
    print("occmax:", occmax)

    Elist = scipy.linspace(Emin, Elim, Elim-Emin+1)
    print("Elist:", Elist)
    print("addTails:", addTails)
    print("saveondb:", saveondb)

    if saveondb:
        if json == False:
            db = database.Database()
        else:
            db = database.Database(dbname="spectraJson.db", useJson=True)

    a = phi4.Phi4()
    a.buildBasis(Emax=Emaxbar, L=L, m=m, k=k, occmax=occmax)
    print("Basis size: ", a.basis[k].size)

    try:
        a.loadMatrix(L=L, Emax=Emaxbar, k=k, occmax=occmax)
        print("matrix loaded")
    except FileNotFoundError:
        print("building matrix...")
        a.buildMatrix(k=k)

    a.setCouplings(0, 0, g)
    # b = finiteVolH.FiniteVolH(a.L, m)
    # g0, g2, g4 = b.directCouplings(g)

    for Emax in Elist:
        print("Emax: ", Emax)

        if addTails:
            cutoff = Emaxbar
        else:
            cutoff = Emax

        if saveondb:
            approxQuery = {"g":g, "L":L, "Emaxbar":cutoff, "Emax":Emax}
            exactQuery = {"k":k, "occmax":occmax}
            if db.getObjList('spec', approxQuery=approxQuery, exactQuery=exactQuery) != []:
                print("Eigenvalues already present")
                continue

        Er = 0
        for ren in ("raw","renlocal"):
            a.renlocal(Emax=cutoff, Er=Er)
            a.computeHamiltonian(Emax=Emax, k=k, ren=ren, addTails=addTails)

            compsize = a.compH.shape[0]
            print("Comp basis size: ", a.compH.shape[0])

            a.computeEigval(k=k, ren=ren, sigma=sigma, neigs=neigs)
            Er = a.vacuumE(ren="raw")

            print("{} vacuum: ".format(ren), a.vacuumE(ren=ren))

            if saveondb:
                if addTails:
                    Emaxbardb = Emaxbar
                else:
                    Emaxbardb = Emax
                db.insert(k=k, Emax=Emax, L=a.L, ren=ren, g=g, spec=a.eigenvalues[ren][k],
                        Emaxbar=cutoff, eigv=a.eigenvectors[ren][k],
                        occmax=occmax, basisSize=compsize, neigs=neigs)


if __name__ == "__main__":
    main(sys.argv)
