import os
import phi4
import sys
import math
import scipy

m = 1.

def main(argv):
    args = "<L> <Emin> <Emax> <k> <?occmax>"
    if len(argv) < 4:
        print("python", argv[0], args)
        sys.exit(-1)

    L = float(argv[1])
    Emin = int(argv[2])
    Emax = int(argv[3])
    k = int(argv[4])
    try:
        occmax = int(argv[5])
    except IndexError:
        occmax = None


    sizelist = []
    nelemlist = []
    nopslist = []

    for E in range(Emin, Emax):
        a = phi4.Phi4()
        a.buildBasis(k=k, Emax=E, m=m, L=L, occmax=occmax)
        print("basis size:", a.basis[k].size)

        a.buildMatrix(k=k)

        nelemlist.append(len(a.V[k][4].M.data))
        nopslist.append(a.nops)
        sizelist.append(a.basis[k].size)
        # print("number of matrix elements:",n)
        # print("sparsity:",n/(a.basis[k].size**2))

    print(list(range(Emin,Emax)))
    print(nelemlist)
    print(nopslist)
    print(sizelist)

if __name__ == "__main__":
    main(sys.argv)
