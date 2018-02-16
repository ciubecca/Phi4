import statefuncs
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


output = "pdf"

plt.style.use('ggplot')
plt.rc('axes', prop_cycle=(cycler('color', ['r', 'g', 'b', 'y'])))
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

params = {'legend.fontsize': 8, 'lines.markersize':2.5, 'lines.marker':"o"}
plt.rcParams.update(params)

plt.rc('axes', prop_cycle=(
    cycler('marker', ['x', 'o', 'v','x'])
    +cycler('linestyle', ['-', '--', ':','-'])
    +cycler('markersize', [4.,2.,2.,2.])
    +cycler('color', ['r','b','g','k'])
    ))


m = 1
# Number of eigenvalues to compute per sector
neigs = 2


argv = sys.argv
if len(argv) < 2:
    print(argv[0], " <L>")
    sys.exit(-1)

L = float(argv[1])

print("L", L)

glist = np.linspace(1, 1.5, 30)
# glist = np.linspace(0.01, 0.1, 30)
print("glist", glist)

xmax = max(glist)+0.03
xmin = 0

ETlist = [10, 15, 20]

vevrawlist = {ET:[] for ET in ETlist}
vevrenlist = {ET:[] for ET in ETlist}
Lambda = {ET:[] for ET in ETlist}

for ET in ETlist:

    a = phi4.Phi4(m, L, 1)
    a.buildBasis(Emax=ET)
    a.computePotential()
    a.setglist(glist=glist)

    b = phi4.Phi4(m, L, -1)
    b.buildBasis(Emax=ET)
    b.computePotential()
    b.setglist(glist=glist)


    a.computeEigval(ET, "raw", neigs=neigs)
    b.computeEigval(ET, "raw", neigs=neigs)
    E0raw = {g: a.eigenvalues[g]["raw"][0] for g in glist}
    E1raw = {g: b.eigenvalues[g]["raw"][0] for g in glist}
    massraw = {g:E1raw[g] - E0raw[g] for g in glist}
    massraw = {g:E1raw[g] - E0raw[g] for g in glist}
    vevrawlist[ET] = [a.vev[g]["raw"] for g in glist]

    a.computeEigval(ET, "renloc", neigs=neigs, eps=E0raw)
    b.computeEigval(ET, "renloc", neigs=neigs, eps=E1raw)
    E0ren = {g: a.eigenvalues[g]["renloc"][0] for g in glist}
    E1ren = {g: b.eigenvalues[g]["renloc"][0] for g in glist}
    massren = {g:E1ren[g] - E0ren[g] for g in glist}
    Lambda[ET] = np.array([E0ren[g] for g in glist])/L
    vevrenlist[ET] = [a.vev[g]["renloc"] for g in glist]


# Estimate of VEV taking derivative of vacuum energy density.
vev2 = dict()
for ET in ETlist:
    vev2[ET] = 2*(Lambda[ET]-glist*np.gradient(Lambda[ET], glist[1]-glist[0]))/(1+3*glist/pi)


plt.figure(1)
for ET in ETlist:
    plt.plot(glist, vevrenlist[ET], label=r"$\langle \phi^2 \rangle$ , Emax={}".format(ET))
    plt.plot(glist, vev2[ET], label=r"$\mathcal{{E}}'$ , Emax={}".format(ET))


plt.figure(1)
# plt.ylim(0,0.12)
plt.xlabel("g")
plt.ylabel(r"$\langle\phi^2\rangle$")
plt.title(r"$L$ = {}".format(L))
plt.legend()
fname = ".{0}".format(output)
s = "VEVvsG_L={}".format(L)
# plt.savefig(s+fname, bbox_inches='tight')
plt.savefig(s+fname)

