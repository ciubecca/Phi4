from statefuncs import Basis

m=1
L=10
Emax=10
k=1

b = Basis.fromScratch(m=m, L=L, Emax=Emax, k=k)
b2 = Basis.fromBasis(b, lambda v: v.energy<=8)
# print(b)
print(b2)
