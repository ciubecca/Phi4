import newphi4
import phi4
import sys
import time
import numpy as np

Emax = 14
L = 10
m = 1

(PotEven, PotOdd) = newphi4.newphi4(ET=Emax,L=L,m=m)

a = phi4.Phi4()

for k in (-1,1):
    a.buildBasis(k=k, Emax=Emax, m=m, L=L)
    print("Lorenzo's basis size :", a.basis[k].size)

    a.buildMatrix(k=k)

print(a.V[1][4].M.nnz)
print(a.V[-1][4].M.nnz)

print("maximal differences of matrices produced by two methods")
print (np.max((PotEven-a.V[1][4].M).todense()), np.max((PotOdd-a.V[-1][4].M).todense()))
