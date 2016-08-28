import os
import phi4
import sys
import math
import scipy

m = 1.

def main(argv):
    args = "<L> <Emax> <k> <?occmax>"
    if len(argv) < 4:
        print("python", argv[0], args)
        sys.exit(-1)

    L = float(argv[1])
    Emax = float(argv[2])
    k = int(argv[3])
    try:
        occmax = int(argv[4])
    except IndexError:
        occmax = None

    a = phi4.Phi4()

    a.buildBasis(k=k, Emax=Emax, m=m, L=L, occmax=occmax)
    print("basis size:", a.basis[k].size)

    a.buildMatrix(k=k)
    n = len(a.V[k][4].M.data)
    print("number of matrix elements:",n)
    print("sparsity:",n/(a.basis[k].size**2))

if __name__ == "__main__":
    main(sys.argv)
