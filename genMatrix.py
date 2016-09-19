import inspect
import os
import phi4
import sys
import math
import scipy

m = 1.

def main(argv):
    args = "<L> <Emax> <?occmax>"
    if len(argv) < 3:
        print("python", argv[0], args)
        sys.exit(-1)

    L = float(argv[1])
    Emax = float(argv[2])
    try:
        occmax = int(argv[3])
    except IndexError:
        occmax = None

    a = phi4.Phi4(m,L,Emax)

    a.buildBasis(Emax=Emax, occmax=occmax)
    print("basis sizes :", a.basis[1].size, a.basis[-1].size)

    # print(a.basis[k])


    a.buildMatrix()
    # a.saveMatrix(k=k)

    # print(a.V[k][4].M.todense())

if __name__ == "__main__":
    main(sys.argv)
