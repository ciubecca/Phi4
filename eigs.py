import statefuncs
import phi4
import renorm
import sys
import scipy
import math
import database

saveondb = False
m = 1.
klist = (1,)
neigs = 10

argv = sys.argv

if len(argv) < 5:
    print(argv[0], " <L> <g> <ET> <EL>")
    sys.exit(-1)

L = float(argv[1])
g = float(argv[2])
ET = float(argv[3])
EL = float(argv[4])

if saveondb:
    db = database.Database()

a = phi4.Phi4(m, L)
a.buildBasis(Emax=ET)


for k in klist:

    print("k=", k)

    if saveondb:
        # Emaxbar == Emax means there are no tails
        approxQuery = {"g":g, "L":L, "ET":ET, "EL":EL}
        exactQuery = {"k":k}
        if db.getObjList('spec', approxQuery=approxQuery, exactQuery=exactQuery) != []:
            print("Eigenvalues already present")
            continue

    print("Basis size: ", a.basis[k].size)

    # try:
        # a.loadMatrix(L=L, Emax=Emaxbar, k=k, occmax=occmax)
        # print("matrix loaded")
    # except FileNotFoundError:
        # print("building matrix...")
        # a.buildMatrix(k=k)

    a.computePotential(k)

    a.setCouplings(0, 0, g)
    # b = finiteVolH.FiniteVolH(a.L, m)
    # g0, g2, g4 = b.directCouplings(g)

    a.computeEigval(k, ET, "raw", neigs=10)
    eps = a.vacuumE("raw")
    print("Raw vacuum:", eps)


    vectorlist = [state for i,state in enumerate(a.basis[k])
            if a.eigenvectors["raw"][1][0][i] > 10**-3]
    basisl = statefuncs.Basis(k, vectorlist, a.basis[k].helper)
    print("Size of selected basis:", basisl.size)

    print("Generating high energy basis...")
    a.genHEBasis(k, basisl, ET, EL)

    print("Computing DeltaH...")
    a.computeHEVs(k, ET, EL)

    a.computeEigval(k, ET, "ren", EL, eps, neigs=10)
    print("Renormalized vacuum:", a.vacuumE("ren"))

    if saveondb:
        # If Emaxbar == Emax it means there are no tails
        db.insert(k=k, ET=ET, EL=EL, L=L, ren="raw", g=g,
                spec=a.eigenvalues["raw"][k], eigv=a.eigenvectors["raw"][k],
                basisSize=a.compSize[k], neigs=neigs)
        db.insert(k=k, ET=ET, EL=EL, L=L, ren="ren", g=g,
                spec=a.eigenvalues["ren"][k], eigv=a.eigenvectors["ren"][k],
                basisSize=a.compSize[k], neigs=neigs, ntails=ntails)

