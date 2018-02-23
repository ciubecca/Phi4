import scipy
import scipy.sparse.linalg
import scipy.sparse

def occn(state):
    """ Computes the occupation number of a state"""
    return sum(Zn for n,Zn in state)

def buildStatePos(basis, helper=None, Erange=None):
    # Dictionary of positions of states
    # Contains also the P-reversed states
    # NOTE: using arrays is much less efficient!

    statePos = {}
    if helper==None:
        helper = basis.helper

    if Erange == None:
        irange = range(basis.size)
    else:
        irange = basis.irange(Erange)

    for i in irange:
        state = basis.repr2List[i]
        statePos[state] = i
        statePos[state[::-1]] = i

    return statePos


class SubmatrixOperator():

    def __init__(self, idxList):
        self.idxList = scipy.array(idxList)

    @classmethod
    def fromSubbasis(self, basis, subbasis):
        statePos = buildStatePos(basis)
        idxList = [statePos[state] for state in subbasis.repr2List]
        return self(idxList)

    @classmethod
    def fromErange(self, basis, Erange):
        return self(basis.irange(Erange))

    def subrows(self, m):
        return m.tocsr()[self.idxList,]

    def subcolumns(self, m):
        return m.tocsc()[:,self.idxList]

    def sub(self, m):
        return self.subrows(self.subcolumns(m))


class MatrixConstructor():
    def __init__(self, basis, lookupbasis, Erange=None):
        """
        basis: basis for the row elements
        lookupbasis: basis for the columns elements
        Emin: minimal energy of the states to be generated
        Emax: maximal energy of the states to be generated
        """
        self.basis = basis
        self.lookupbasis = lookupbasis

        if basis.helper.nmax > lookupbasis.helper.nmax:
            self.helper = basis.helper
        else:
            self.helper = lookupbasis.helper

        if Erange==None:
            Erange = (lookupbasis.Emin, lookupbasis.Emax)
        self.Erange = Erange

        self.statePos = buildStatePos(lookupbasis, self.helper, Erange)

    def buildMatrix(self, Vlist, ignKeyErr=False, idxList=None, sumTranspose=False):
        """
        Vlist: list of oscillators
        ignKeyErr: whether LookupError when generating a state should be ignored (to be used
        when the lookupbasis does not contain all the states between Emin and Emax)
        idxList: list of the row indices to be computed
        sumTranspose: whether the tranpose of the matrix should be added and the diagonal
        subtracted
        """

        basis = self.basis
        lookupbasis = self.lookupbasis
        Erange = self.Erange
        statePos = self.statePos
        helper = self.helper

        if idxList==None:
            idxList = range(basis.size)

        # Will construct the sparse matrix in the COO format and then convert it to CSC
        data = []
        row = []
        col = []

        for V in Vlist:
            for i in idxList:

                # XXX Check added to avoid generating states with higher occupation number than cutoff
                if lookupbasis.occmax!=None and V.nc-V.nd+occn(basis.stateList[i]) > lookupbasis.occmax:
                    continue

                colpart, datapart = \
                    V.computeMatrixElements(basis,i,lookupbasis, Erange=Erange,
                            statePos=statePos, helper=helper, ignKeyErr=ignKeyErr)
                data += datapart
                col += colpart
                row += [i]*len(colpart)

        # XXX Does this sum duplicate entries?
        V = scipy.sparse.coo_matrix((data,(row,col)), shape=(basis.size,lookupbasis.size))

        if sumTranspose:
            # Add the matrix to its transpose and subtract the diagonal
            diag_V = scipy.sparse.spdiags(V.diagonal(),0,basis.size,basis.size).tocsc()
            return (V+V.transpose()-diag_V).tocsc()
        else:
            return V.tocsc()

    def __del__(self):
        del self.statePos
