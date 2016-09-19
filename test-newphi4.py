from newphi4Lorenzo import newPhi4
import sys

a = newPhi4()

args = "<L> <Emax>"
if len(sys.argv) < 3:
    print sys.argv[0], args
    sys.exit(-1)

L = float(sys.argv[1])
ET = float(sys.argv[2])
m = 1

a.buildBasisAndMatrix(L=L, ET=ET, m=m)

print "basis size", len(a.basis[1])

print a.potential[1].todense()
print a.potential[-1].todense()
