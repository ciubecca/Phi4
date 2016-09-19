import inspect
import os
import phi4
import sys
import math
import scipy
from statefuncs import *

m = 1.

ET = 10

def main(argv):
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

    a = phi4.Phi4(m,L,Ebar)

    nmax =  a.info.nmax
    print("nmax", nmax)

    a.buildBasis(Emax=Emax, occmax=occmax)

    # print(a.basis[1])
    print("basis sizes :", a.basis[1].size, a.basis[-1].size)



    # print(a.basis[k])

    subbasis = a.basis[1].sub(lambda v: occn(v)<=4)
    a.computeDH2(1, subbasis, Emax, Ebar)


    print(a.basisH[1].energyList)
    print("HE basis size", a.basisH[1].size)
    print("Emax", a.basisH[1].Emin)
    print("Emin", a.basisH[1].Emax)


    print("Estimated memory use of basis:", 16*(2*nmax+1)*a.basisH[1].size/(10^6), "M")

    # a.buildMatrix()
    # a.saveMatrix(k=k)

    # print(a.V[k][4].M.todense())

if __name__ == "__main__":
    main(sys.argv)
