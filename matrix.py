from profile_support import *
import scipy
import scipy.sparse.linalg
import scipy.sparse

def buildStatePos(basis):
    # Dictionary of positions of states
    # TODO Contains also the P-reversed states
    # NOTE: using arrays is much less efficient!

    statePos = {}
    helper = basis.helper
    irange = range(basis.size)

    for i in irange:
        state = basis.repr2List[i]
        statePos[state] = i

        # XXX To implement
        # statePos[state[::-1]] = i

    return statePos



class MatrixConstructor():
    def __init__(self, basis):
        """
        basis: basis for the row and column elements
        """
        self.basis = basis
        self.statePos = buildStatePos(basis)

    @profile
    def buildMatrix(self, Vlist, ignKeyErr=False, sumTranspose=True):
        """
        Vlist: list of oscillators
        ignKeyErr: whether LookupError when generating a state should be ignored (to be used
        when the lookupbasis does not contain all the states between Emin and Emax)
        idxList: list of the row indices to be computed
        sumTranspose: whether the tranpose of the matrix should be added and the diagonal
        subtracted
        """

        basis = self.basis
        statePos = self.statePos
        helper = self.basis.helper

        idxList = range(basis.size)

        # Will construct the sparse matrix in the COO format and then convert it to CSC
        data = []
        row = []
        col = []

        for V in Vlist:
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
