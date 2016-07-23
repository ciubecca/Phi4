import finiteVolH
import phi4
import renorm
import sys
import scipy
import math
import database

json = False

def main(argv):
    if len(argv) < 2:
        print(argv[0], " <L> <g>")
        return -1

    L = float(argv[1])
    g = float(argv[2])

    # Hardcoded parameters
    m = 1.
    Emaxbar = 30
    Elist = scipy.linspace(6, 29, 24)
    occmax = 4
    sigma = -30.
    neigs = 1
    k = 1
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

        # Emaxbar == Emax means there are no tails
        approxQuery = {"g":g, "L":L, "Emaxbar":Emax, "Emax":Emax}
        exactQuery = {"k":k}
        if db.getObjList('spec', approxQuery=approxQuery, exactQuery=exactQuery) != []:
            print("Eigenvalues already present")
            continue

        Er = 0
        for ren in ("raw","renlocal"):
            # This is different wrt with tails
            a.renlocal(Emax=Emax, Er=Er)
            a.computeHamiltonian(Emax=Emax, k=k, ren=ren, Er=Er, addTails=False)

            compsize = a.compH.shape[0]
            print("Comp basis size: ", a.compH.shape[0])

            a.computeEigval(k=k, ren=ren, sigma=sigma, neigs=neigs)
            Er = a.vacuumE(ren="raw")

            print("{} vacuum: ".format(ren), a.vacuumE(ren=ren))

            # If Emaxbar == Emax it means there are no tails
            db.insert(k=k, Emax=Emax, Emaxbar=Emax, L=a.L, ren=ren, g=g, spec=a.eigenvalues[ren][k],
                    eigv=a.eigenvectors[ren][k], occmax=occmax, basisSize=compsize, neigs=neigs)


if __name__ == "__main__":
    main(sys.argv)
