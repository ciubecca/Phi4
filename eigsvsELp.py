import statefuncs
import phi4
import renorm
import sys
import scipy
import math
import database

# Suppose we want to generate all the eigenvalues with L=10, g=1 with truncation energies
# ET = 10, 10.5, 11, ..., 20
# Then we should call this file as:
# eigsvsE.py 10 1 10 20

# Whether we should save the results in the database data/spectra.db
saveondb = True
# saveondb = False
m = 1
# List of parity quantum numbers
klist = (1,-1)
# Minimum overlap with the raw vacuum for selecting a state in the tails
minoverlap = 10**(-2)
# Ratio between ELpp and ELp
ratio3 = 1.5
# Ratio between EL and ET
ratio2 = 2
neigs = 10

# List of all the contributions to DH3. Sequentially, we add DH3<<, DH3<> and DH3>>
tlist = ((False,False,False),(True,True,False),(True,True,True))



argv = sys.argv
if len(argv) < 6:
    print(argv[0], "<L> <g> <ET> <ELpmin> <ELpmax>")
    sys.exit(-1)

L = float(argv[1])
g = float(argv[2])
ET = float(argv[3])
ELpmin = float(argv[4])
ELpmax = float(argv[5])

# EL is chosen as the maximum ELpp
EL = ratio2*ET

ELplist = scipy.linspace(ELpmin, ELpmax, (ELpmax-ELpmin)*2+1)
print("ELplist:", ELplist)

print("minoverlap:", minoverlap)
print("ELpp/ELp:", ratio3)
print("EL/ET:", ratio2)

if saveondb:
    db = database.Database()

a = phi4.Phi4(m, L)
a.buildBasis(Emax=ET)


for k in klist:

    # Compute the potential matrices in the low-energy space below ET
    a.computePotential(k)

    print("k=", k)
    print("Full basis size: ", a.basis[k].size)

    a.setCouplings(g4=g)
    print("g=", g)

    print("Computing raw eigenvalues for highest cutoff")
    a.computeEigval(k, ET, "raw")

    # Select a set of tails and construct a Basis object
    vectorlist = [state for i,state in enumerate(a.basis[k])
        if abs(a.eigenvectors["raw"][k][0][i]) > minoverlap]
    basisl = statefuncs.Basis(k, vectorlist, a.basis[k].helper)
    print("Total number of tails:", basisl.size)


    print("Generating high energy basis...")
    # Generate the high-energy "selected" basis by passing a set of tails
    # and a maximum cutoff EL
    a.genHEBases(k, basisl, EL=EL, ELpp=ratio3*ELpmax)
    print("Size of HE basis:", a.basisH[k].size)

    a.computeLEVs(k)


    print("Computing high energy matrices...")
# Compute the matrices VLH, VHL, VHH, for the highest local cutoff ELmax.
# Later we will be varying EL, therefore taking submatrices of these.
# Computing VHH is expensive
    a.computeHEVs(k)



# Compute the raw eigenvalues for cutoff ET
    a.computeEigval(k, ET, "raw")
    print("Raw vacuum:", a.eigenvalues["raw"][k][0])
    eps = a.eigenvalues["raw"][k][0]


# Compute "local" renormalized eigenvalues for cutoff ET
# Since we are passing EL=ET to the method call, the matrices VHL, VHH will be computed
# only in the local approximation
    a.computeEigval(k, ET, "renloc", EL=ET, eps=eps)
    print("Local ren vacuum:", a.eigenvalues["renloc"][k][0])
    eps = a.eigenvalues["renloc"][k][0]



    if True in (t[2] for t in tlist):
        a.calcVV3(ELplist, eps)

    for t in tlist:

        nonloc3mix, loc3mix, loc3 = t

        for ELp in ELplist:

            ELpp = ratio3*ELp
            print("ELp={}, ELpp={}".format(ELp,ELpp))

            a.computeEigval(k, ET, "rentails", EL=EL, ELp=ELp, ELpp=ELpp, eps=eps,
                    loc2=True,loc3=loc3, loc3mix=loc3mix, nonloc3mix=nonloc3mix, neigs=neigs)
            print("Non-Local ren vacuum:", a.eigenvalues["rentails"][k][0])

            print("Number of tails:", a.ntails)


            if saveondb:
                datadict = dict(k=k, ET=ET, L=L, ren="rentails", g=g, minoverlap=minoverlap,
                    EL=EL, ELp=ELp, ELpp=ELpp, ntails=a.ntails, eps=eps, neigs=neigs,
                    loc3=loc3, loc3mix=loc3mix, nonloc3mix=nonloc3mix, basisSize=a.compSize)

                db.insert(datadict=datadict, spec=a.eigenvalues["rentails"][k])

    del a.VLH[k]
