import statefuncs
import phi4
import renorm
import sys
import scipy
import math
import database


# Whether we should save the results in the database data/spectra.db
saveondb = True
# saveondb = False
m = 1
# List of parity quantum numbers
klist = (1,-1)
# Maximum number of tails
# NOTE The actual number of tails can increase for higher ET
maxntails = 300

# Ratio between EL and ET
ratioELET = 3
# Ratio between ELp and ET
ratioELpET = 1.5
# Ratio between ELpp and ELp
ratioELppELp = 1.5

neigs = 10

argv = sys.argv
if len(argv) < 5:
    print(argv[0], "<L> <g> <ETmin> <ETmax>")
    sys.exit(-1)

L = float(argv[1])
g = float(argv[2])
ETmin = float(argv[3])
ETmax = float(argv[4])

ELmax = ratioELET*ETmax
ELpmax = ratioELpET*ETmax
ELppmax = ratioELppELp*ELpmax

print("EL/ET:", ratioELET)
print("ELp/ET:", ratioELpET)
print("ELpp/ELp", ratioELppELp)


ETlist = scipy.linspace(ETmin, ETmax, (ETmax-ETmin)*2+1)
print("ETlist:", ETlist)

if saveondb:
    db = database.Database()

a = phi4.Phi4(m, L)
a.buildBasis(Emax=ETmax)

for k in klist:
    print("k=", k)

    print("Full basis size: ", a.basis[k].size)
    a.computePotential(k)

# Computing the high energy bases dimension for highest Emax
    a.setCouplings(g4=g)
    a.computeEigval(k, ETmax, "raw")


    # Select a set of tails and construct a Basis object, ordered in overlap with
    # the vacuum
    vectorlist = [state for i,state in sorted(enumerate(a.basis[k]), key=lambda x:
            -abs(a.eigenvectors["raw"][k][0][x[0]]))][:maxntails]
    basisl = statefuncs.Basis(k, vectorlist, a.basis[k].helper)
    print("ntails:", basisl.size)


    print("Generating high energy basis for highest Emax...")
    # Generate the high-energy "selected" basis by passing a set of tails
    # and a maximum cutoff EL
    a.genHEBases(k, basisl, EL=ELmax, ELp=ELpmax, ELpp=ELppmax)
    print("Size of HE basis", a.basisH[k].size)

    a.computeLEVs(k)

    print("Computing high energy matrices...")
# Compute the matrices VLH, VHL, VHH, for the highest local cutoff ELmax.
# Later we will be varying EL, therefore taking submatrices of these.
# Computing VHH is expensive
    a.computeHEVs(k)


    for ET in ETlist:

        print("ET", ET)

        EL = ratioELET*ET
        ELp = ratioELpET*ET
        ELpp = ratioELppELp*ELp

        a.computeEigval(k, ET, "raw", neigs=neigs)
        print("Raw vacuum:", a.eigenvalues["raw"][k][0])
        eps = a.eigenvalues["raw"][k][0]


        a.computeEigval(k, ET, "renloc", eps=eps, neigs=neigs)
        print("Local ren vacuum:", a.eigenvalues["renloc"][k][0])
        eps = a.eigenvalues["renloc"][k][0]

        a.calcVV3([ELp], eps)


        a.computeEigval(k, ET, "rentails", EL=EL, ELp=ELp, ELpp=ELpp, eps=eps, neigs=neigs)
        print("Non-Local ren vacuum:", a.eigenvalues["rentails"][k][0])

        print("Number of tails:", a.ntails)


        if saveondb:
            datadict = dict(k=k, ET=ET, L=L, ren="raw", g=g, neigs=neigs,
                    basisSize=a.compSize)
            db.insert(datadict=datadict, spec=a.eigenvalues["raw"][k])


            datadict = dict(k=k, ET=ET, L=L, ren="renloc", g=g, eps=eps, neigs=neigs,
                    basisSize=a.compSize)
            db.insert(datadict=datadict, spec=a.eigenvalues["renloc"][k])


            datadict = dict(k=k, ET=ET, L=L, ren="rentails", g=g, EL=EL, ELp=ELp, ELpp=ELpp,
                    ntails=a.ntails, eps=eps, neigs=neigs, basisSize=a.compSize,
                    tailsComputedAtET=ETmax, maxntails=maxntails)
            db.insert(datadict=datadict, spec=a.eigenvalues["rentails"][k])


    del a.VLH[k]
