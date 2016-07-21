from statefuncs import Basis
from phi4 import Phi4

m = 1
L = 10
Emax = 12
k = 1

b = Basis.fromScratch(m=m, L=L, Emax=Emax, k=k, occmax=4)

b = Basis.fromBasis(b, lambda v: v.energy <= 8)
print(b)
# print(b2)

# a = Phi4()
# a.buildBasis(Emax = Emax, L=L, m=m, k=k)
