import finiteVolH
import phi4
import renorm
import sys
import scipy
import math
import database

json = False
occmax = 3
m = 1.
sigma = -30.
neigs = 1
klist = (-1,)

def main(argv):
    if len(argv) < 4:
        print(argv[0], " <L> <Emax> <g>")
        return -1

    L = float(argv[1])
    Emax = float(argv[2])
    g = float(argv[3])

    if json == False:
        db = database.Database()
    else:
        db = database.Database(dbname="data/spectraJson.db", useJson=True)

    a = phi4.Phi4()

    for k in klist:
        print("k=", k)
        a.buildBasis(Emax=Emax, L=L, m=m, k=k, occmax=occmax)
        print("Basis size: ", a.basis[k].size)

        a.buildMatrix(k=k)
        a.setCouplings(0, 0, g)
        # b = finiteVolH.FiniteVolH(a.L, m)
        # g0, g2, g4 = b.directCouplings(g)

        # Emaxbar == Emax means there are no tails
        approxQuery = {"g":g, "L":L, "Emaxbar":Emax, "Emax":Emax}
        exactQuery = {"k":k, "occmax":occmax}
        if db.getObjList('spec', approxQuery=approxQuery, exactQuery=exactQuery) != []:
            print("Eigenvalues already present")
            continue

        Er = 0
        for ren in ("raw","renlocal"):
            # This is different wrt with tails
            a.renlocal(Emax=Emax, Er=Er)
            a.computeHamiltonian(Emax=Emax, k=k, ren=ren, Er=Er, addTails=False)

            compsize = a.compH.shape[0]
            if ren == "raw": print("Comp basis size: ", a.compH.shape[0])

            a.computeEigval(k=k, ren=ren, sigma=sigma, neigs=neigs)
            # Er = a.vacuumE(ren="raw")

            if k==1: print("{} vacuum: ".format(ren), a.vacuumE(ren=ren))

            # If Emaxbar == Emax it means there are no tails
            db.insert(k=k, Emax=Emax, Emaxbar=Emax, L=a.L, ren=ren, g=g,
                    spec=a.eigenvalues[ren][k], eigv=a.eigenvectors[ren][k],
                    occmax=occmax, basisSize=compsize, neigs=neigs)

if __name__ == "__main__":
    main(sys.argv)
