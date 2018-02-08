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

output = "pdf"

plt.style.use('ggplot')
plt.rc('axes', prop_cycle=(cycler('color', ['r', 'g', 'b', 'y'])))
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

params = {'legend.fontsize': 8, 'lines.markersize':2.5, 'lines.marker':"o"}
plt.rcParams.update(params)

plt.rc('axes', prop_cycle=(
    cycler('marker', ['x', 'o', 'v'])
    +cycler('linestyle', ['-', '--', ':'])
    +cycler('markersize', [5.,3.,3.])
    ))

m = 1
# Number of eigenvalues to compute per sector
neigs = 2


argv = sys.argv
if len(argv) < 3:
    print(argv[0], " <L> <ET>")
    sys.exit(-1)

L = float(argv[1])
ET = float(argv[2])


print("L, ET", L, ET)

glist = np.linspace(0.1, 2, 20)
print("glist", glist)

xmax = max(glist)+0.03
xmin = 0


def main():

    a = phi4.Phi4(m, L, 1)
    a.buildBasis(Emax=ET)
    a.computePotential()
    a.setglist(glist=glist)

    b = phi4.Phi4(m, L, -1)
    b.buildBasis(Emax=ET)
    b.computePotential()
    b.setglist(glist=glist)

    a.computeLEVs()
    b.computeLEVs()


    a.computeEigval(ET, "raw", neigs=neigs)
    b.computeEigval(ET, "raw", neigs=neigs)
    E0raw = {g: a.eigenvalues[g]["raw"][0] for g in glist}
    E1raw = {g: b.eigenvalues[g]["raw"][0] for g in glist}
    massraw = {g:E1raw[g] - E0raw[g] for g in glist}
    # vevraw = np.array([a.vev[g]["raw"] for g in glist])
    print("Raw vacuum:", E0raw)
    print("Raw mass", massraw)
    # print("raw vev:", vevraw)

    a.computeEigval(ET, "renloc", neigs=neigs, eps=E0raw)
    b.computeEigval(ET, "renloc", neigs=neigs, eps=E1raw)
    E0ren = {g: a.eigenvalues[g]["renloc"][0] for g in glist}
    E1ren = {g: b.eigenvalues[g]["renloc"][0] for g in glist}
    massren = {g:E1ren[g] - E0ren[g] for g in glist}
    # vevraw = np.array([a.vev[g]["raw"] for g in glist])
    print("ren vacuum:", E0ren)
    print("ren mass", massren)
    # print("raw vev:", vevraw)

    # plt.figure(1, figsize=(4, 2.5), dpi=300)
    plt.figure(1)

    massrawlist = np.array([massraw[g] for g in glist])
    massrenlist = np.array([massren[g] for g in glist])
    plt.plot(glist, massrawlist, label="ET raw")
    plt.plot(glist, massrenlist, label="ET ren")
    plt.xlim(xmin, xmax)
    # plt.ylim(0,1)

    plt.title(r"$E_T$ = {} , $L$ = {}".format(ET, L))
    plt.legend()
    fname = ".{0}".format(output)
    s = "LC_ET={}_L={}".format(ET,L)
    # plt.savefig(s+fname, bbox_inches='tight')
    plt.savefig(s+fname)


main()
