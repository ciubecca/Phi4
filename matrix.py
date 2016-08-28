import scipy
import scipy.sparse.linalg
import scipy.sparse
import scipy.interpolate

# TODO use metaclasses?
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

    def sub(self, subBasisI, subBasisJ):
        """ This extracts a submatrix given a subspace of
        the initial vector space, both for rows and columns
        """
        rows = [self.basisI.lookup(state)[1]  for state in subBasisI]
        columns = [self.basisJ.lookup(state)[1]  for state in subBasisJ]
        return Matrix(subBasisI, subBasisJ,
                self.M.tocsr()[scipy.array(rows),].
                tocsc()[:,scipy.array(columns)])

    def transpose(self):
        return Matrix(self.basisJ, self.basisI, self.M.transpose())
    def inverse(self):
        return Matrix(self.basisJ, self.basisI, scipy.sparse.linalg.inv(self.M.tocsc()))
