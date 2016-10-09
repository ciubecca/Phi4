import statefuncs
import phi4
import renorm
import sys
import scipy
import math
import database

saveondb = True
# saveondb = False
m = 1
neigs = 10
klist = (1,)
minoverlap = 10**(-2)

argv = sys.argv
if len(argv) < 6:
    print(argv[0], "<L> <g> <ETmin> <ETmax> <EL-ET>")
    sys.exit(-1)

L = float(argv[1])
g = float(argv[2])
ETmin = float(argv[3])
ETmax = float(argv[4])
ELETdiff = float(argv[5])

ETlist = scipy.linspace(ETmin, ETmax, ETmax-ETmin+1)
print("ETlist:", ETlist)

if saveondb:
    db = database.Database()

a = phi4.Phi4(m, L)
a.buildBasis(Emax=ETmax)

a.setCouplings(g4=g)

for k in klist:

    a.computePotential(k)

    print("k=", k)
    print("Full basis size: ", a.basis[k].size)


    print("Computing raw eigenvalues for highest cutoff")
    a.computeEigval(k, ETmax, "raw", neigs=neigs)
    eps = a.vacuumE("raw")

    print(a.eigenvectors["raw"][1][0])

    vectorlist = [state for i,state in enumerate(a.basis[k])
        if abs(a.eigenvectors["raw"][1][0][i]) > minoverlap]
    basisl = statefuncs.Basis(k, vectorlist, a.basis[k].helper)
    print("Total number of tails:", basisl.size)

    ELmax = ETmax + ELETdiff

    print("Generating high energy basis...")
    a.genHEBasis(k, basisl, ELmax)
    print("Size of HE basis:", a.basisH[k].size)

    print("Computing high energy matrices...")
    a.computeHEVs(k, ELmax)


    for i,ET in enumerate(ETlist):
        EL = ET + ELETdiff
        print("ET={}, EL={}".format(ET,EL))

        if saveondb:
            # Emaxbar == Emax means there are no tails
            approxQuery = {"g":g, "L":L, "ET":ET, "EL":EL}
            exactQuery = {"k":k, "ren":"ren"}
            if db.getObjList('spec', approxQuery=approxQuery,
                    exactQuery=exactQuery) != []:
                print("Eigenvalues already present")
                continue

        a.computeEigval(k, ET, "raw", neigs=neigs)
        print("Raw vacuum:", a.vacuumE("raw"))

        a.computeEigval(k, ET, "ren", EL, eps, neigs=neigs)
        print("Renormalized vacuum:", a.vacuumE("ren"))

        print("Number of tails:", a.ntails)

        if saveondb:
            # If Emaxbar == Emax it means there are no tails
            db.insert(k=k, ET=ET, L=L, ren="raw", g=g,
                    spec=a.eigenvalues["raw"][k], eigv=a.eigenvectors["raw"][k],
                    basisSize=a.compSize, EL=EL, neigs=neigs)
            db.insert(k=k, ET=ET, L=L, ren="ren", g=g,
                    spec=a.eigenvalues["ren"][k], eigv=a.eigenvectors["ren"][k],
                    basisSize=a.compSize, neigs=neigs, EL=EL, ntails=a.ntails)
