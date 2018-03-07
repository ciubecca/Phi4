import statefuncs
from scipy import stats
import phi4
import matplotlib.pyplot as plt
from matplotlib import rc
from cycler import cycler
import renorm
import sys
import scipy
import math
import numpy as np
from numpy import pi, sqrt, log, exp
import gc


output = "pdf"

plt.style.use('ggplot')
plt.rc('axes', prop_cycle=(cycler('color', ['r', 'g', 'b', 'y'])))
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

params = {'legend.fontsize': 8, 'lines.markersize':2.5, 'lines.marker':"o"}
plt.rcParams.update(params)

plt.rc('axes', prop_cycle=(
    cycler('marker', ['x', 'o', 'v','x','o'])
    +cycler('linestyle', ['-', '--', ':','-','--'])
    +cycler('markersize', [4.,2.,2.,2.,2.])
    +cycler('color', ['r','b','g','k','c'])
    ))


m = 1
# Number of eigenvalues to compute per sector
neigs = 5
neigs = 2

argv = sys.argv
if len(argv) < 2:
    print(argv[0], " <g>")
    sys.exit(-1)

g = float(argv[1])

glist = [g]
print("g", g)


Llist = np.linspace(5,12,8)
print("Llist", glist)

ETlist = [15, 18, 20]
# ETlist = [10, 15]
# ETlist = [10]

for ET in ETlist:

    phi01 = []
    Gap = []

    for L in Llist:

        print("ET: ", ET)

        a = phi4.Phi4(m, L, 1)
        b = phi4.Phi4(m, L, -1)
        a.buildBasis(Emax=ET)
        b.buildBasis(Emax=ET)
        a.computePotential(other=b.basis)

        # print(a.V[1])

        b.computePotential()
        a.setglist(glist=glist)
        b.setglist(glist=glist)

        print("Basis size even:", len(a.basis.stateList))
        print("Basis size odd:", len(b.basis.stateList))


        a.computeEigval(ET, "raw", neigs=neigs)
        b.computeEigval(ET, "raw", neigs=neigs)
        E0raw = {g: a.eigenvalues[g]["raw"][0] for g in glist}
        E1raw = {g: b.eigenvalues[g]["raw"][0] for g in glist}
        massraw = E1raw[glist[0]] - E0raw[glist[0]]

        a.computeEigval(ET, "renloc", neigs=neigs, eps=E0raw)
        b.computeEigval(ET, "renloc", neigs=neigs, eps=E1raw)
        E0ren = {g: a.eigenvalues[g]["renloc"][0] for g in glist}
        E1ren = {g: b.eigenvalues[g]["renloc"][0] for g in glist}
        massren = E1ren[glist[0]] - E0ren[glist[0]]

        Gap.append(massren)

        phi01.append(2/L*massren*np.inner(a.eigenvectors[glist[0]]["renloc"][0], a.V[1].dot(b.eigenvectors[glist[0]]["renloc"][0]))**2)

    Gap = np.array(Gap)
    phi01 = np.array(phi01)

    plt.figure(1)
    plt.plot(Llist, phi01, label=r"$\langle 0 \mid \phi \mid 1 \rangle$ , Emax={}".format(ET))
    # plt.figure(2)
    # plt.plot(Llist/(Gap[ET]**2), phi01[ET], label=r"Emax={}".format(ET))



plt.figure(1)
plt.xlim(min(Llist),max(Llist))
plt.xlabel("L")
plt.ylabel(r"$\langle 0 \mid \phi \mid 1 \rangle^2$")
s = "phi01vsL_g={}".format(g)
plt.title(r"$g$ = {}".format(g))
plt.legend()
fname = ".{0}".format(output)
# plt.savefig(s+fname, bbox_inches='tight')
plt.savefig(s+fname)
