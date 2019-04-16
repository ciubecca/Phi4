# Classes and functions related to matrix construction


from profile_support import *
import scipy
import scipy.sparse.linalg
import scipy.sparse

tol = 10**-10

def checkZero(M):
    try:
        assert abs(M).max() < tol
    except AssertionError as e:
        print(abs(M).max())
        raise e
    return

def checkSymmetric(M):
    checkZero(M-M.transpose())
    return

def submatrix(V, subidx):
    """ Return the submatrix for states within smaller cutoffs """
    return (V.tocsc()[:,subidx]).tocsr()[subidx,]

def subrows(V, subidx):
    return V.tocsr()[subidx,]

def subcolumns(V, subidx):
    return V.tocsc()[:,subidx]


def buildMatrix(basis, Vlist, destbasis=None, idxList=None,
        ignKeyErr=False, sumTranspose=True):
    """
    Vlist: list of oscillators
    ignKeyErr: whether LookupError when generating a state should be ignored (to be used
    when the lookupbasis does not contain all the states between Emin and Emax)
    idxList: list of the row indices to be computed
    sumTranspose: whether the tranpose of the matrix should be added and the diagonal
    subtracted
    """

    if destbasis==None:
        destbasis = basis

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
