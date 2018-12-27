from statefuncs import *
from symmetry import *
from profile_support import *
import me
from me import *
from oscillators import *
from matrix import *


def genVHl(basis, subidx, basisH, L):

        Vlist = V4OpsSelectedFull(basis, basisH.helper, subidx)

        return buildMatrix(basis, Vlist, idxList=subidx,
                destbasis=basisH, sumTranspose=False)*L**2


def genHEBases(bases, tailidx, EL, ELp):
        """ Generate a high-energy basis from a set of tails
        basis: Basis containing the set of tail states
        EL: maximal energy of the generated basis for DH2
        ELp: maximal energy of the generated basis for DH3. Only these states
        will be stored in representation 1
        """

# Usually EL > ELp
        Emax = max(EL, ELp)

        helper = Helper(L=bases[1].helper.L, Emax=Emax,
                    Lambda=bases[1].helper.Lambda, m=bases[1].helper.m)


        ret = {}
        for k in (-1,1):
            basis = bases[k]
            subidx = tailidx[k]

            # Generate all the operators between the selected states and the states
            # in the range [0, Emax]
            # XXX Assuming that V2 does not generate different states
            Vlist = V4OpsSelectedFull(basis, helper, subidx)

            statePos = {}
            stateList = []
            ncomp = []

            i = 0
            for V in Vlist:
                for v in V.yieldBasis(basis, subidx, Emax):
                    if bytes(tuple(v)) not in statePos:

                        va = np.array(v)
                        transStates = {bytes(tuple(m.dot(va))) for m in helper.transfMat}

                        stateList.append(bytes(tuple(v)))
                        ncomp.append(len(transStates))

                        for s in transStates:
                            statePos[s] = i
                        i += 1

# Basis of selected states with energy <= Emax. We only need to save
# states in the type 1 representation (most memory consuming) for states
# with energy <= ELp, or ET when ELp=None
            ret[k] = Basis(k, stateList, helper, statePos, ncomp, repr1=False,
                    repr1Emax=max(ELp or 0, basis.Emax))

        return ret



def V4OpsSelectedFull(basis, helper, idxList=None):
    """ Selected set of oscillators of the full V4 operator between some selected states
    basis: basis which is acted upon
    idxList: subset of indices of states which are acted upon
    helper: Helper object of the destination basis
    """

    oscEnergy = helper.oscEnergy
    Emax = helper.Emax

    if idxList == None:
        idxList = range(basis.size)

    opsList = []
    allowedWnPairs = None

    for nd in (0,1,2,3,4):
        nc = 4-nd

        dlists = gendlistsfromBasis(basis, idxList, helper, nd, 4)

        if nd==2:
            totpairsmomenta = set((k1[0]+k2[0],k1[1]+k2[1]) for k1,k2 in dlists)
            allowedWnPairs = helper.genMomentaPairs(totpairsmomenta)

        oscList = []
        for dlist in dlists:
            clists = [clist for clist in
                        createClistsV4(helper, dlist, nc, allowedWnPairs) if
                        oscEnergy(clist) <= Emax+tol]
            oscList.append((dlist, clists))

        opsList.append(LocOperator(oscList,nd,nc,helper=helper))

    return opsList


def gendlistsfromBasis(basis, idxList, helper, nd, ntot):
    ret = set()

    for i in idxList:
        state = basis.stateList[i]
        ret.update(gendlists(state=state, nd=nd, ntot=ntot, helper=helper))

    return ret


def createClistsV4(helper, dlist, nc, allowedWnPairs=None):

    if len(dlist) != 4-nc:
        raise ValueError
    clists = []

    if nc==0:
        clists.append(())

    elif nc==1:
        clists.append((sum(map(np.array, dlist)),))

    elif nc==2:
        k1,k2 = dlist
        k12 = (k1[0]+k2[0],k1[1]+k2[1])

        for k3,k4 in allowedWnPairs[k12]:
            clists.append((k3,k4))

    elif nc==3:
        (k1,) = dlist
        clists = helper.genMomenta3sets(k1)

    elif nc==4:
        clists = []
        clists =  helper.genMomenta4sets()

    return clists
