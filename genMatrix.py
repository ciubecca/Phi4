import inspect
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

    a = phi4.Phi4(m,L,Emax)

    a.buildBasis(k=k, Emax=Emax, occmax=occmax)
    print("basis size :", a.basis[k].size)

    # print(a.basis[k])


    a.buildMatrix(k=k)
    # a.saveMatrix(k=k)

    # print(a.V[k][4].M.todense())

if __name__ == "__main__":
    main(sys.argv)
