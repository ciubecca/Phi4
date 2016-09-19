from newphi4Lorenzo import newPhi4
import sys

a = newPhi4()

args = "<L> <Emax> <k>"
if len(sys.argv) < 4:
    print sys.argv[0], args
    sys.exit(-1)

L = float(sys.argv[1])
ET = float(sys.argv[2])
k = int(sys.argv[3])
m = 1

a.buildBasisAndMatrix(L=L, ET=ET, m=m)

# print "basis size", len(a.basis[k])

# print a.potential[k].todense()
# print a.potential[-1].todense()
