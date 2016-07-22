from statefuncs import Basis
from phi4 import Phi4

m = 1
L = 10
Emaxbar = 24
Emax = 22
k = 1
g = 0.1
ren = "renlocal"

# b = Basis.fromScratch(m=m, L=L, Emax=Emax, k=k, occmax=4)
# b = Basis.fromBasis(b, lambda v: v.energy <= 8)
# print(b)

a = Phi4()
a.buildBasis(Emax = Emaxbar, L=L, m=m, k=k, occmax=4)
print("Basis size: ", a.basis[1].size)

a.buildMatrix(k=k)
a.setCouplings(0, 0, g)
print("Matrices built.")

a.renlocal(Emax=Emax, Er=0.)
# a.computeEigval(k=1, Emax = 10, ren="renlocal")
# print(a.eigenvalues["renlocal"][1])

a.computeHamiltonian(k=k, Emax=Emax, ren=ren)
print("Computational basis size: ", a.Hfull.shape[0])
a.computeEigval(ren,k)

# a.computeHamiltonian(k=1, Emax=Emax, ren="renlocal")
# a.computeEigval(k=1, Emax=Emax, ren="renlocal")

print(a.eigenvalues[ren][k])
