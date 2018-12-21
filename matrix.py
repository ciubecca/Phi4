from profile_support import *
import scipy
import scipy.sparse.linalg
import scipy.sparse

def submatrix(V, subidx):
    """ Return the submatrix for states within smaller cutoffs """
    return (V.tocsc()[:,subidx]).tocsr()[subidx,]

def subrows(V, subidx):
    return V.tocsr()[subidx,]

def subcolumns(V, subidx):
    return V.tocsc()[:,subidx]



class MatrixConstructor():
    def __init__(self, basis, destbasis=None):
        """
        basis: basis for the row elements
        destbasis: destination basis
        """
        self.basis = basis
        if destbasis==None:
            self.destbasis = basis
        else:
            self.destbasis = destbasis


    def buildMatrix(self, Vlist, idxList=None, ignKeyErr=False, sumTranspose=True):
        """
        Vlist: list of oscillators
        ignKeyErr: whether LookupError when generating a state should be ignored (to be used
        when the lookupbasis does not contain all the states between Emin and Emax)
        idxList: list of the row indices to be computed
        sumTranspose: whether the tranpose of the matrix should be added and the diagonal
        subtracted
        """

        basis = self.basis
        destbasis = self.destbasis

        if idxList==None:
            idxList = range(basis.size)

        # Will construct the sparse matrix in the COO format and then convert it to CSC
        data = []
        row = []
        col = []

        for V in Vlist:
            for i,idx in enumerate(idxList):
                colpart, datapart = \
                    V.computeMatrixElements(basis, idx, destbasis, ignKeyErr=ignKeyErr)
                    # statePos=statePos, ignKeyErr=ignKeyErr)
                data += datapart
                col += colpart
                row += [i]*len(colpart)

        # XXX Does this sum duplicate entries?
        # XXX Check
        V = scipy.sparse.coo_matrix((data,(row,col)),
                shape=(len(idxList), destbasis.size))

        if sumTranspose:
            # Add the matrix to its transpose and subtract the diagonal
            diag_V = scipy.sparse.spdiags(V.diagonal(),0,basis.size,
                    basis.size).tocsc()
            return (V+V.transpose()-diag_V).tocsc()
        else:
            return V.tocsc()
