from phi4 import *
from sys import argv, exit
from math import factorial
from statefuncs import *
import copy
from scipy import sparse

sym = True
m = 1
neigs = 4

if len(argv) < 5:
    print("{} <L> <Emax> <Lambda> <g4> [<g2>]".format(argv[0]))
    exit(1)

L = float(argv[1])
Emax = float(argv[2])
Lambda = float(argv[3])
g4 = float(argv[4])
try:
    g2 = float(argv[5])
except IndexError:
    g2 = 0

# print("gap for g4=0 :".format(sqrt(1+2*g2)))


print("Computing basis...")
bases = Basis.fromScratch(m, L, Emax, Lambda)

for k in (-1,1):
    print("k={}, L={}, Emax={}, Lambda={}, size={}".format(k, L, Emax, Lambda, len(bases[k])))

eigs = {}

debug = False

def printMatrix(m):
    md = m.todense()
    for i in range(m.shape[0]):
        for j in range(i+1, m.shape[0]):
            md[i,j] = 0
    print(sparse.csr_matrix(md))

Vlist = None
V22 = None

for k in (-1,1):
    print("Computing k={} matrices...".format(k))
    a = Phi4(bases[k])
    Vlist, V22 = a.computePotential(Vlist, V22)

    a.setg(0, g2, g4/(factorial(4)))

    print("Diagonalizing matrix...".format(k))
    a.computeEigval(neigs=neigs)
    eigs[k] = a.eigval

    # if k==-1 and debug:
        # # print("state[15]: ", toCanonical(a.basis.stateList[15]))
        # # print("state[3]: ", toCanonical(a.basis.stateList[3]))
        # # print(a.basis.energyList)

        # for i,s in enumerate(a.basis.stateList):
            # if i==0 or i==6:
                # print("{}: {}".format(i,s))
                # print("e: {}".format(a.basis.energyList[i]))

        # assert abs((a.V[4]-a.V[4]).transpose()).max() < 10**-5
        # printMatrix(a.V[4])
        # np.savetxt("matrix.csv", a.V[4].todense().reshape(1,a.basis.size**2), delimiter=',')
        # pass


print("L={}, Emax={}, Lambda={}, g2={}, g4={}".format(L, Emax, Lambda, g2, g4))
print("Vacuum: ", eigs[1][0])
print("Even spectrum: ", eigs[1][1:]-eigs[1][0])
print("Odd  spectrum: ", eigs[-1]-eigs[1][0])
