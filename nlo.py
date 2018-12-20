from statefuncs import *
from symmetry import *

def genHEBasis(self, basis, subidx, EL, ELp):
        """ Generate a high-energy basis from a set of tails
        k: parity quantum number
        basis: Basis containing the set of tail states
        EL: maximal energy of the generated basis for DH2
        ELp: maximal energy of the generated basis for DH3. Only these states
        will be stored in representation 1
        """

        m = basis.helper.m
        Lambda = basis.helper.Lambda

# Usually EL > ELp
        Emax = max(EL, ELp)
# Helper function of the new basis
        helper = Helper(Emax=Emax, Lambda=Lambda, m=m)

        # Generate all the operators between the selected states and the states
        # in the range [0, Emax]
        # XXX Assuming that V2 does not generate different states
        Vlist = V4OpsSelectedFull(basis, helper, subidx)

        statePos = {}
        stateList = []

        i = 0
        for V in Vlist:
            for v in V.yieldBasis(basis, subidx, Emax):
                if v not in statePos:
                    stateList.append(v)
                    for m in helper.transfMat
                        vt = tuple(np.dot(m, v))
                        statePos[vt] = i
                    i += 1

# Basis of selected states with energy <= Emax. We only need to save
# states in the type 1 representation (most memory consuming) for states
# with energy <= ELp, or ET when ELp=None
        self.basisH = Basis(self.k, stateList, helper, statePos, repr1=False,
                repr1Emax=max(ELp or 0, basis.Emax))
