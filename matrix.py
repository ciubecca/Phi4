import scipy
import scipy.sparse.linalg
import scipy.sparse

class Matrix():
    """ Matrix with specified state bases for row and column indexes.
    This class is useful to easily extract submatrices """
    def __init__(self, basisI, basisJ, M):
        self.basisI = basisI
        self.basisJ = basisJ
        self.M = M
        self.check()

    def check(self):
        if self.M.shape != (self.basisI.size, self.basisJ.size):
            raise ValueError('Matrix shape inconsistent with given bases')

    def __add__(self, other):
        """ Sum of matrices """
        return Matrix(self.basisI, self.basisJ, self.M+other.M)
    def __sub__(self, other):
        return Matrix(self.basisI, self.basisJ, self.M-other.M)

    def __rmul__(self, other):
        """ Multiplication of matrix with matrix or number"""
        if(other.__class__ == self.__class__):
            return Matrix(other.basisI, self.basisJ, other.M*self.M)
        else:
            return Matrix(self.basisI, self.basisJ, self.M*float(other))

    def __mul__(self, other):
        """ Multiplication of matrix with matrix or number"""
        if(other.__class__ == self.__class__):
            return Matrix(self.basisI, other.basisJ, self.M*other.M)
        else:
            return Matrix(self.basisI, self.basisJ, self.M*float(other))

    def to(self, form):
        """ Format conversion """
        return Matrix(self.basisI, self.basisJ, self.M.asformat(form))

    def subIndex(self, Ibounds, Jbounds):
        """ Returns a submatrix bz bounding row and columns in a given range.
        It should be faster than the more general sub() routine """
        raise RuntimeError("Not yet implemented")
        return

    def sub(self, subBasisI, subBasisJ):
        """ This extracts a submatrix given a subspace of
        the initial vector space, both for rows and columns
        """

        if subBasisI == self.basisI:
            columns = [self.basisJ.lookup(state) for state in subBasisJ]
            return Matrix(subBasisI, subBasisJ, self.M.tocsc()[:,scipy.array(columns)])

        elif subBasisJ == self.basisJ:
            rows = [self.basisI.lookup(state) for state in subBasisI]
            return Matrix(subBasisI, subBasisJ, self.M.tocsr()[scipy.array(rows),])

        else:
            rows = [self.basisI.lookup(state) for state in subBasisI]
            columns = [self.basisJ.lookup(state) for state in subBasisJ]

            return Matrix(subBasisI, subBasisJ, self.M.tocsr()[scipy.array(rows),].
                                                tocsc()[:,scipy.array(columns)])

    def transpose(self):
        return Matrix(self.basisJ, self.basisI, self.M.transpose())
    def inverse(self):
        return Matrix(self.basisJ, self.basisI, scipy.sparse.linalg.inv(self.M.tocsc()))
