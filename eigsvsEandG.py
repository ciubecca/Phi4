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
neigs = 1
klist = (-1,1)
minoverlap = 10**(-2)

ETmin = 10
ETmax = 20
gmin = 0.5
gmax = 3.0

ETlist = scipy.linspace(ETmin, ETmax, (ETmax-ETmin)*2+1)
print("ETlist:", ETlist)
glist = scipy.linspace(gmin, gmax, (gmax-gmin)*2+1)
print("glist:", glist)


def el(ET):
    return ET*2.5

argv = sys.argv
if len(argv) < 2:
    # print(argv[0], "<L> <g> <ETmin> <ETmax> <EL-ET>")
    print(argv[0], "<L>")
    sys.exit(-1)

L = float(argv[1])
# ELETdiff = float(argv[5])

if saveondb:
    db = database.Database()

a = phi4.Phi4(m, L)
a.buildBasis(Emax=ETmax)


for k in klist:

    a.computePotential(k)

    print("k=", k)
    print("Full basis size: ", a.basis[k].size)


    print("Computing raw eigenvalues for highest cutoff")
    a.computeEigval(k, ETmax, "raw", neigs=neigs)


    vectorlist = [state for i,state in enumerate(a.basis[k])
        if abs(a.eigenvectors["raw"][k][0][i]) > minoverlap]
    basisl = statefuncs.Basis(k, vectorlist, a.basis[k].helper)
    print("Total number of tails:", basisl.size)

    ELmax = el(ETmax)

    print("Generating high energy basis...")
    a.genHEBasis(k, basisl, ELmax)
    print("Size of HE basis:", a.basisH[k].size)

    print("Computing high energy matrices...")
    a.computeHEVs(k, ELmax)


    for ET in ETlist:

        for g in glist:

            a.setCouplings(g4=g)
            print("g=", g)

            EL = el(ET)
            # EL = ELmax
            print("ET={}, EL={}".format(ET,EL))

            if saveondb:
                # Emaxbar == Emax means there are no tails
                approxQuery = {"g":g, "L":L, "ET":ET, "EL":EL}
                exactQuery = {"k":k, "ren":"rentails"}
                if db.getObjList('spec', approxQuery=approxQuery,
                        exactQuery=exactQuery) != []:
                    print("Eigenvalues already present")
                    continue

            a.computeEigval(k, ET, "raw", neigs=neigs)
            print("Raw vacuum:", a.eigenvalues["raw"][k][0])
            eps = a.eigenvalues["raw"][k][0]

            if saveondb:
                # If Emaxbar == Emax it means there are no tails
                db.insert(k=k, ET=ET, L=L, ren="raw", g=g,
                        spec=a.eigenvalues["raw"][k], eigv=a.eigenvectors["raw"][k],
                        basisSize=a.compSize, neigs=neigs)


            a.computeEigval(k, ET, "renloc", ET, eps, neigs=neigs)
            print("Local ren vacuum:", a.eigenvalues["renloc"][k][0])
            eps = a.eigenvalues["renloc"][k][0]

            if saveondb:
                db.insert(k=k, ET=ET, L=L, ren="renloc", g=g,
                        spec=a.eigenvalues["renloc"][k], eigv=a.eigenvectors["renloc"][k],
                        basisSize=a.compSize, neigs=neigs, EL=ET, eps=eps)


            a.computeEigval(k, ET, "rentails", EL, eps, neigs=neigs)
            print("Non-Local ren vacuum:", a.eigenvalues["rentails"][k][0])

            print("Number of tails:", a.ntails)

            if saveondb:
                db.insert(k=k, ET=ET, L=L, ren="rentails", g=g,
                        spec=a.eigenvalues["rentails"][k], eigv=a.eigenvectors["rentails"][k],
                        basisSize=a.compSize, neigs=neigs, EL=EL, ntails=a.ntails, eps=eps)
