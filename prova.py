from statefuncs import Basis
from phi4 import Phi4

m = 1
L = 10
Emax = 8
k = 1

# b = Basis.fromScratch(m=m, L=L, Emax=Emax, k=k, occmax=4)
# b = Basis.fromBasis(b, lambda v: v.energy <= 8)
# print(b)

a = Phi4()
a.buildBasis(Emax = Emax, L=L, m=m, k=k)
a.buildMatrix(k=1)
a.setCouplings(0, 0, 0.1)
a.renlocal(Emax=10, Er=-0.1)
# a.computeEigval(k=1, Emax = 10, ren="renlocal")
# print(a.eigenvalues["renlocal"][1])


a.computeHamiltonian(k=1, Emax=6, ren="raw")
