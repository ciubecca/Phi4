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
neigs = 20
neigs = 2


argv = sys.argv
if len(argv) < 2:
    print(argv[0], " <L>")
    sys.exit(-1)

L = float(argv[1])

print("L", L)

# glist = np.linspace(1, 1.5, 30)
glist = np.linspace(0.1, 2, 20)
# glist = np.linspace(0, 0.1, 2)
print("glist", glist)

xmax = max(glist)+0.03
xmin = 0

ETlist = [15, 20, 22]
# ETlist = [10]

vevrawlist = {ET:[] for ET in ETlist}
vevrenlist = {ET:[] for ET in ETlist}
Lambda = {ET:[] for ET in ETlist}
phi01 = {ET:[] for ET in ETlist}
commutator = {ET:[] for ET in ETlist}
intSpec = {ET:[] for ET in ETlist}

Gap = {}

for ET in ETlist:

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
    massraw = {g:E1raw[g] - E0raw[g] for g in glist}

    a.computeEigval(ET, "renloc", neigs=neigs, eps=E0raw)
    b.computeEigval(ET, "renloc", neigs=neigs, eps=E1raw)
    E0ren = {g: a.eigenvalues[g]["renloc"][0] for g in glist}
    E1ren = {g: b.eigenvalues[g]["renloc"][0] for g in glist}
    massren = {g:E1ren[g] - E0ren[g] for g in glist}
    Lambda[ET] = np.array([E0ren[g] for g in glist])/L

    # print([(a.V[1]/a.L).dot(b.eigenvectors[g]["renloc"][0]) for g in glist])

    # print(massren[glist[10]])
    # print(b.eigenvectors[glist[10]]["renloc"][0])
    # print(b.basis.stateList[0])

    Gap[ET] = np.array([massren[g] for g in glist])

    phi01[ET] = 2/L*np.array([massren[g]*np.inner(a.eigenvectors[g]["renloc"][0], a.V[1].dot(b.eigenvectors[g]["renloc"][0]))**2 for g in glist])

    v1 = {g: a.eigenvectors[g]["renloc"][0] for g in glist}
    v2 = {g: b.eigenvectors[g]["renloc"][0] for g in glist}

    commutator[ET] = 2/L**2*np.array([massren[g]*np.inner(v1[g], a.V[1].dot(v2[g]))*np.inner(v1[g], a.P[1].dot(v2[g]))
              for g in glist])

    for g in glist:
        x = 0
        for e, v in zip(b.eigenvalues[g]["renloc"], b.eigenvectors[g]["renloc"]):
            gap = e - E0ren[g]
            if gap > 5*massren[g] or gap < 3*massren[g]:
                continue
            x +=  2*np.inner(a.eigenvectors[g]["renloc"][0],(a.V[1]).dot(v))**2

        intSpec[ET].append(x)

    gc.collect()


for ET in ETlist:
    plt.figure(1)
    plt.plot(glist, phi01[ET], label=r"$\langle 0 \mid \phi \mid 1 \rangle$ , Emax={}".format(ET))
    plt.figure(2)
    plt.plot(glist, intSpec[ET], label=r"$I(5 \mu) - I(3 \mu)$, Emax={}".format(ET))
    plt.figure(3)
    plt.plot(glist, commutator[ET], label=r"Emax={}".format(ET))
    plt.figure(4)
    plt.plot(glist/(Gap[ET]**2), phi01[ET], label=r"Emax={}".format(ET))

# plt.figure(4)
# plt.show()


plt.figure(1)
plt.xlim(0,max(glist))
plt.xlabel("g")
plt.ylabel(r"$\langle 0 \mid \phi \mid 1 \rangle$")
s = "phi01vsG_L={}".format(L)
plt.title(r"$L$ = {}".format(L))
plt.legend()
fname = ".{0}".format(output)
# plt.savefig(s+fname, bbox_inches='tight')
plt.savefig(s+fname)

plt.figure(2)
plt.xlim(0,max(glist))
plt.xlabel("g")
plt.ylabel(r"$I(5 \mu) - I(3 \mu)$")
s = "SpecDensityvsG_L={}".format(L)
plt.title(r"$L$ = {}".format(L))
plt.legend()
fname = ".{0}".format(output)
# plt.savefig(s+fname, bbox_inches='tight')
plt.savefig(s+fname)


plt.figure(3)
plt.xlim(0,max(glist))
plt.xlabel("g")
plt.ylabel(r"$\langle 0 | [\phi, \pi ] | 0 \rangle$")
s = "CommutatorvsG_L={}".format(L)
plt.title(r"$L$ = {}".format(L))
plt.legend()
fname = ".{0}".format(output)
# plt.savefig(s+fname, bbox_inches='tight')
plt.savefig(s+fname)


plt.figure(4)
plt.xlim(0,max(glist))
plt.xlabel("g")
plt.ylabel(r"$\langle 0 | [\phi, \pi ] | 0 \rangle$")
s = "Matt_L={}".format(L)
plt.title(r"$L$ = {}".format(L))
plt.legend()
fname = ".{0}".format(output)
# plt.savefig(s+fname, bbox_inches='tight')
plt.savefig(s+fname)
