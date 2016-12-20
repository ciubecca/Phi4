import scipy
import scipy.sparse.linalg
import scipy.sparse


def buildStatePos(basis, helper=None):
    # Dictionary of positions of states
    # Contains also the P-reversed states
    # NOTE: using arrays is much less efficient!
    statePos = {}
    if helper==None:
        helper = basis.helper

    for i,state in enumerate(basis.stateList):
        statePos[tuple(helper.torepr2(state))] = i
        statePos[tuple(helper.torepr2(state)[::-1])] = i

    return statePos


class SubmatrixOperator():

    def __init__(self, idxList):
        self.idxList = scipy.array(idxList)

    def fromSubbasis(self, basis, subbasis):
        statePos = buildStatePos(basis)
        idxList = [statePos[state] for state in subbasis]
        return self(idxList)

    def fromEmax(self, basis, Emax):
        return self(range(basis.imax(Emax)))

    def subrows(self, m):
        return m.tocsr()[self.idxList,]

    def subcolumns(self, m):
        return m.tocsc()[:,self.idxList]

    def sub(self, m):
        return self.subrows(self.subcolumns(m))
