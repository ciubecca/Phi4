import finiteVolH
import phi4
import renorm
import sys
import scipy
import math
import database

json = False

# Hardcoded parameters
m = 1.
Emaxbar = 20
Elist = scipy.linspace(5, 19, 15)
occmax = 6
sigma = -30.
neigs = 1
k = 1

def main(argv):
    if len(argv) < 3:
        print(argv[0], " <L> <g>")
        return -1

    L = float(argv[1])
    g = float(argv[2])

    print(Elist)

    if json == False:
        db = database.Database()
    else:
        db = database.Database(dbname="spectraJson.db", useJson=True)

    a = phi4.Phi4()
    a.buildBasis(Emax=Emaxbar, L=L, m=m, k=k, occmax=occmax)
    print("Basis size: ", a.basis[k].size)

    a.buildMatrix(k=k)
    a.setCouplings(0, 0, g)
    # b = finiteVolH.FiniteVolH(a.L, m)
    # g0, g2, g4 = b.directCouplings(g)

    for Emax in Elist:
        print("Emax: ", Emax)

        approxQuery = {"g":g, "L":L, "Emaxbar":Emaxbar, "Emax":Emax}
        exactQuery = {"k":k}
        if db.getObjList('spec', approxQuery=approxQuery, exactQuery=exactQuery) != []:
            print("Eigenvalues already present")
            continue

        Er = 0
        for ren in ("raw","renlocal"):
            a.renlocal(Emax=Emaxbar, Er=Er)
            a.computeHamiltonian(Emax=Emax, k=k, ren=ren, addTails=True, Er=Er)

            compsize = a.compH.shape[0]
            print("Comp basis size: ", a.compH.shape[0])

            a.computeEigval(k=k, ren=ren, sigma=sigma, neigs=neigs)
            Er = a.vacuumE(ren="raw")

            print("{} vacuum: ".format(ren), a.vacuumE(ren=ren))

            db.insert(k=k, Emax=Emax, L=a.L, ren=ren, g=g, spec=a.eigenvalues[ren][k],
                    Emaxbar=Emaxbar, eigv=a.eigenvectors[ren][k],
                    occmax=occmax, basisSize=compsize, neigs=neigs)


if __name__ == "__main__":
    main(sys.argv)
