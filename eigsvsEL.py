import statefuncs
import phi4
import renorm
import sys
import scipy
import math
import database

dbname = "spectravsEL.db"
saveondb = True
test = True

# Ratio between ELp and ET
ratioELpET = 2
# Ratio between ELpp and ELp
ratioELppELp = 1.5

m = 1
k = 1
neigs = 6

db = database.Database(dbname)

argv = sys.argv
if len(argv) < 6:
    print(argv[0], "<L> <g> <ET> <ELmin> <ELmax>")
    sys.exit(-1)

L = float(argv[1])
g = float(argv[2])
ET = float(argv[3])
ELmin = float(argv[4])
ELmax = float(argv[5])

ELp = ratioELpET*ET
ELpp = ratioELpELpp*ELp
print("ELp:", ELp, "ELpp", ELpp)

ELlist = scipy.linspace(ELmin, ELmax, (ELmax-ELmin)*2+1)
print("ELlist:", ELlist)



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

a.genHEBasis(EL=ELmax, ELp=ELp, ELpp=ELpp)
print("Size of HE basis:", a.basisH.size)

a.computeHEVs()

a.computeEigval(ET, "renloc", eps=epsraw, neigs=neigs)
eps = {g: a.eigenvalues[g]["renloc"][0] for g in glist}
print("Local ren vacuum:", eps)

a.calcVV3(ELp, eps, test=test)

for EL in ELlist:
    print("EL={}".format(EL))

    for loc2 in (True, False):
        a.computeEigval(ET, "rentails", EL=EL, ELp=ELp, ELpp=ELpp, eps=eps,
            neigs=neigs, loc2=loc2, loc3=True)

        if saveondb:
            datadict = dict(k=k, ET=ET, L=L, ren="rentails", g=g, EL=EL, ELp=ELp,
                    ELpp=ELpp, ntails=a.ntails, eps=eps[g], neigs=neigs,
                    basisSize=a.compSize, finiteL=True, loc2=loc2)
            db.insert(datadict=datadict, spec=a.eigenvalues[g]["rentails"])
