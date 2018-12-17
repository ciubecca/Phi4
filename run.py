import random
from phi4 import *
from statefuncs import *

Emax= 7
L= 5
occmax = 4
m = 1

basesFull = Basis.fromScratch(m, L, Emax)
bases = Basis.fromScratch(m, L, Emax, occmax=occmax)

k = -1
sub = [s for s in basesFull[k].stateList if occn(s) <= occmax]

print("sub size: {}".format(len(sub)))
print("statelist size: {}".format(bases[k].size))

e = bases[k].

for i in range(len(sub)):
    if sub[i] != bases[k].stateList[i]:
        print(sub[i])
        print(helper.
        print(bases[k].stateList[i])
        print(sub[i+1])
        exit(1)
# assert bases[k].size == len(sub)


