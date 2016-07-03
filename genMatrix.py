import inspect
import os
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
os.sys.path.insert(0,parentdir)

import phi1234
import sys
import math
import scipy

def main(argv):
    
    if len(argv) < 3:
        print "python genMatrix.py <L> <Emax>"
        sys.exit(-1)
    
    L=float(argv[1])
    Emax=float(argv[2])
    
    m = 1.

    a = phi1234.Phi1234()

    a.buildFullBasis(k=1,Emax=Emax,m=m,L=L)
    a.buildFullBasis(k=-1,Emax=Emax,m=m,L=L)

    print "K=1 basis size :", a.fullBasis[1].size
    print "K=-1 basis size :", a.fullBasis[-1].size

    fstr = "Emax="+str(a.fullBasis[1].Emax)+"_L="+str(a.L)

    a.buildMatrix()
    a.saveMatrix(fstr)

if __name__ == "__main__":
    main(sys.argv)
