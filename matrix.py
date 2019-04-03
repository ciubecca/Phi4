from profile_support import *
import scipy
import scipy.sparse.linalg
import scipy.sparse


class MatrixConstructor():
    """ Class used to construct a matrix (or part of it) from a list of oscillators """

    def __init__(self, basis):
        """
        basis: basis for the row and column elements
        """
        self.basis = basis
        self.statePos = basis.statePos

    def buildMatrix(self, Vlist, ignKeyErr=False, sumTranspose=True):
        """
        Vlist: list of oscillators
        ignKeyErr: whether to ignore LookupError when generating a state which is not in the column basis (to be used
        when the column basis does not contain all the states between Emin and Emax)
        idxList: list of the row indices to be computed
        sumTranspose: whether to add to the final result the transpose minus the diagonal
        """

        basis = self.basis
        statePos = self.statePos
        helper = self.basis.helper

        idxList = range(basis.size)

        # Will construct the sparse matrix in the COO format and then convert it to CSC
        data = []
        row = []
        col = []

        # Cicle over list of oscillators
        for V in Vlist:
            # Cycle over row indices
            for i in idxList:
                colpart, datapart = \
                    V.computeMatrixElements(basis, i, statePos=statePos, ignKeyErr=ignKeyErr)
                data += datapart
                col += colpart
                row += [i]*len(colpart)

        # XXX Does this sum duplicate entries?
        V = scipy.sparse.coo_matrix((data,(row,col)), shape=(basis.size,basis.size))

        if sumTranspose:
            # Add the matrix to its transpose and subtract the diagonal
            diag_V = scipy.sparse.spdiags(V.diagonal(),0,basis.size,basis.size).tocsc()
            return (V+V.transpose()-diag_V).tocsc()
        else:
            return V.tocsc()

    def __del__(self):
        del self.statePos
