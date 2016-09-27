import inspect
import os
import phi4
import sys
import math
import scipy
from statefuncs import *

k = 1
m = 1.

argv = sys.argv

args = "<L> <Emax> <Ebar> <?occmax>"
if len(argv) < 4:
    print("python", argv[0], args)
    sys.exit(-1)

L = float(argv[1])
Emax = float(argv[2])
Ebar = float(argv[3])
try:
    occmax = int(argv[4])
except IndexError:
    occmax = None

a = phi4.Phi4(m,L)


a.buildBasis(Emax=Emax, occmax=occmax)

nmax =  a.helper.nmax
# print("nmax", nmax)


# print(a.basis[k])
# print("basis size :", a.basis[k].size)

# print(a.basis[1][-1])
# print(a.basis[1].stateDlists[2][-1])

a.buildMatrix()
# print(a.basis[k])

# print(a.basis[k])


vset = [
[],
[(0, 2)],
[(-1, 1), (1, 1)],
[(-1, 1), (0, 2), (1, 1)],
# [(-1, 2), (0, 1), (1, 2)],
[(-2, 1), (-1, 1), (1, 1), (2, 1)]
]

# vset = [
# [],
# [(0,2)],
# [(-1,1),(1,1)]
# ]

subbasis = Basis(k, vset, a.helper, gendlist=True)
print("subbasis size:", subbasis.size)

a.computeDH2(k, subbasis, Emax, Ebar)


# print(a.basisH[1].energyList)
print("HE basis size", a.basisH[k].size)
# print("Emax", a.basisH[1].Emin)
# print("Emin", a.basisH[1].Emax)

# for n in range(2,9,2):
    # print("n=", n)
    # print("States:")
    # print(a.basisH[k].sub(lambda x: occn(x)==n).size)

# for i,v in enumerate(a.basisH[k].stateList):
    # print(a.basisH[k].energyList[i])
    # print(v)

# for state in a.basisH[k].stateList[:10]:
    # print(state)
# print(a.basisH[k].energyList[:10])

# print("Estimated memory use of basis:", 16*(2*nmax+1)*a.basisH[1].size/(10^6), "M")

# a.saveMatrix(k=k)

# print(a.V[k][4].M.todense())

# print(a.V[k][4].M.todense())


b = phi4.Phi4(m,L)
b.buildBasis(Emax=Ebar, occmax=occmax)
b.buildMatrix()

subIndexList = [i for i in range(b.basis[k].size) if b.basis[k].stateList[i] in vset]

V = b.V[1][4].M.todok()

HEindexSet = set()
for i in subIndexList:
    for j,v in enumerate(b.basis[k].stateList):
        if b.basis[k].energyList[j] > Emax and V[i,j] != 0:
            HEindexSet.add(j)

# print(len(HEindexSet))

# print(len([b.basis[k].lookup(state) for state in a.basisH[k]]))

# for state in a.basisH[k]:
    # print(state)

HEindexSet2 = set()
for i, state in enumerate(a.basisH[k]):
    try:
        HEindexSet2.add(b.basis[k].lookup(state))
    except KeyError:
        print("KeyError:", state)
        # print(a.basisH[k].energyList[i])



# HEindexList2 = []
# d = {}
# for j, state in enumerate(a.basisH[k]):
    # i = b.basis[k].lookup(state)
    # if i in HEindexList2:
        # print("Already present:", d[i], j)
        # print("Already present:", a.basisH[k].stateList[d[i]], b.basis[k].stateList[i], state)
    # d[i] = j
    # HEindexList2.append(i)


print(len(HEindexSet))
print(len(HEindexSet2))


# diff = HEindexSet2 - HEindexSet

# print(diff)
