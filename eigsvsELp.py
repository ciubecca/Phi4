import statefuncs
import phi4
import renorm
import sys
import scipy
import math
import database

dbname = "data/spectravsEL.db"
saveondb = True
test = False

# Ratio between ELp and ET
ratioELET = 3
# Ratio between ELpp and ELp
ratioELppELp = 1.5

m = 1
k = 1
neigs = 6

db = database.Database(dbname)

argv = sys.argv
if len(argv) < 6:
    print(argv[0], "<L> <g> <ET> <ELpmin> <ELpmax>")
    sys.exit(-1)

L = float(argv[1])
g = float(argv[2])
ET = float(argv[3])
ELpmin = float(argv[4])
ELpmax = float(argv[5])

EL = ratioELET*ET
ELppmax = ratioELppELp*ELpmax
print("EL:", EL)

ELplist = scipy.linspace(ELpmin, ELpmax, (ELpmax-ELpmin)*2+1)
print("ELplist:", ELplist)

a = phi4.Phi4(m, L, k)
a.buildBasis(Emax=ET)
print("basis size: ", a.basis.size)

a.computePotential()
glist = [g]
a.setglist([g])

a.computeEigval(ET, "raw", neigs=neigs)
epsraw = {g: a.eigenvalues[g]["raw"][0] for g in glist}

basisl = a.basis
a.computeLEVs(basisl, loc3=True)

a.genHEBasis(EL=EL, ELp=ELpmax, ELpp=ELppmax)
print("Size of HE basis:", a.basisH.size)

a.computeHEVs()

a.computeEigval(ET, "renloc", eps=epsraw, neigs=neigs)
eps = {g: a.eigenvalues[g]["renloc"][0] for g in glist}
print("Local ren vacuum:", eps)


for ELp in ELplist:
    ELpp = ratioELppELp*ELp
    print("ELp={}, ELpp={}".format(ELp,ELpp))

    a.calcVV3(ELp, eps, test=test)

    for (nonloc3mix, loc3mix, loc3) in ((False, False, False),(True,True,True)):
        a.computeEigval(ET, "rentails", EL=EL, ELp=ELp, ELpp=ELpp, eps=eps,
            neigs=neigs, nonloc3mix=nonloc3mix, loc3mix=loc3mix, loc3=loc3,
            memsave=False)

        if saveondb:
            datadict = dict(k=k, ET=ET, L=L, ren="rentails", g=g, EL=EL, ELp=ELp,
                    ELpp=ELpp, ntails=a.ntails, eps=eps[g], neigs=neigs,
                    basisSize=a.compSize, finiteL=True, test=test,
                    loc2=True, nonloc3mix=nonloc3mix, loc3mix=loc3mix, loc3=loc3)
            db.insert(datadict=datadict, spec=a.eigenvalues[g]["rentails"])
