import statefuncs
from sklearn.linear_model import LinearRegression
import phi4
import matplotlib.pyplot as plt
from matplotlib import rc
from cycler import cycler
import renorm
import sys
import scipy
import math
import numpy as np
from numpy import log, sqrt, exp
from scipy import stats

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
neigs = 2

argv = sys.argv

L = 5
ET = 10
g = 1

print("L, ET, g", L, ET, g)

a = phi4.Phi4(m, L, 1)
b = phi4.Phi4(m, L, -1)
a.buildBasis(Emax=ET)
b.buildBasis(Emax=ET)
a.computePotential(other=b.basis)
b.computePotential()
a.setglist(glist=[g])
b.setglist(glist=[g])

print(a.stateList[0])
print(a.stateList[1])

phi01[ET] = np.array([np.inner(a.eigenvectors[g]["raw"][0],(a.V[1]/a.L).dot(b.eigenvectors[g]["raw"][0]))**2
        for g in glist])

    print(a.basis.stateList[0])
    print(a.basis.stateList[1])

    # This should not be zero
    print(a.V[2].todense()[0,1])

    b = phi4.Phi4(m, L, -1)
    b.buildBasis(Emax=ET)
    b.computePotential()
    b.setglist(glist=glist)

    # a.computeLEVs()
    # b.computeLEVs()


    a.computeEigval(ET, "raw", neigs=neigs)
    b.computeEigval(ET, "raw", neigs=neigs)
    E0raw = {g: a.eigenvalues[g]["raw"][0] for g in glist}
    E1raw = {g: b.eigenvalues[g]["raw"][0] for g in glist}
    massraw = {g:E1raw[g] - E0raw[g] for g in glist}
    vevraw = {g: a.vev[g]["raw"] for g in glist}
    # print("Raw vacuum:", E0raw)
    # print("Raw mass", massraw)
    # print("raw vev:", vevraw)

    a.computeEigval(ET, "renloc", neigs=neigs, eps=E0raw)
    b.computeEigval(ET, "renloc", neigs=neigs, eps=E1raw)
    E0ren = {g: a.eigenvalues[g]["renloc"][0] for g in glist}
    E1ren = {g: b.eigenvalues[g]["renloc"][0] for g in glist}
    massren = {g:E1ren[g] - E0ren[g] for g in glist}
    vevren = {g: a.vev[g]["renloc"] for g in glist}
    # print("ren vacuum:", E0ren)
    # print("ren mass", massren)
    # print("ren vev:", vevren)



main()
