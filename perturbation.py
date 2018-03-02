from statefuncs import *
from oscillators import *
from matrix import *
from scipy import stats
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
neigs = 2

argv = sys.argv
if len(argv) < 2:
    print(argv[0], " <L>")
    sys.exit(-1)

L = float(argv[1])

ETlist = scipy.linspace(10,60,10)

print("L", L)

# Perturbative vacuum state
basis0 = Basis.fromScratch(m=m, L=L, k=1, Emax=max(ETlist), occmax=0)
# Basis of states with 4 particles
basis4 = Basis.fromScratch(m=m, L=L, k=1, Emax=max(ETlist), occmax=4, occmin=4)

# Generate partial 8-particle basis states by acting on all 4-particle states
print("Generating 8 particle basis")
Vlist = V4OpsSelectedFull(basis4, max(ETlist), ndmax=3)
vectorset = set()

for V in Vlist:
    for v in V.yieldBasis(basis4, max(ETlist)):
# Do not add the vacuum
        if sum(v)==0:
            print("Error")
            continue
        # Don't add twice states connected by parity inversion
        if v not in vectorset and v[::-1] not in vectorset:
            vectorset.add(v)


helper = Vlist[0].helper
basis8 = Basis(1, vectorset, helper, repr1=False, repr1Emax=0)
print("Basis size:", basis8.size)

c1 = MatrixConstructor(basis0, basis4)
Vlist = V4OpsSelectedFull(basis0, max(ETlist))
V04 = c1.buildMatrix(Vlist, sumTranspose=False)*L

c2 = MatrixConstructor(basis4, basis4)
Vlist = V4OpsSelectedHalf(basis4, max(ETlist))
V44 = c2.buildMatrix(Vlist, sumTranspose=True)*L

print("Constructing the matrix")
c3 = MatrixConstructor(basis4, basis8)
Vlist = V4OpsSelectedFull(basis4, max(ETlist))
V48 = c3.buildMatrix(Vlist, sumTranspose=False)*L


Delta2 = []
Delta3 = []
Delta4 = []

print("Computing corrections")

###### SECOND ORDER #############
for ET in ETlist:
    prop4 = basis4.propagator(0,0,ET)
    M = V04*prop4*V04.transpose()/L
    Delta2.append(M.todense()[0,0])

###### THIRD ORDER #############
    M = V04*prop4*V44*prop4*V04.transpose()/L
    Delta3.append(M.todense()[0,0])

###### FOURTH ORDER #############
    prop8 = basis8.propagator(0,0,ET)
    M = (V04*prop4*V48*prop8*V48.transpose()*prop4*V04.transpose()-L*Delta2[-1]*V04*(prop4**2)*V04.transpose())/L
    Delta4.append(M.todense()[0,0])


print(Delta2)
print(Delta3)
print(Delta4)

xlist = log(ETlist)/(ETlist**2)
plt.plot(xlist, Delta2, label="L={}".format(L))
# plt.plot(ETlist, Delta2, label="L={}".format(L))
plt.xlim(0, max(xlist))

# plt.xlabel(r"$1/E_T^2$")
plt.xlabel(r"$\log E_T/E_T^2$")
plt.ylabel(r"$\Lambda$")
plt.savefig("E0pert2.pdf")
plt.clf()

xlist = 1/(ETlist**3)
plt.plot(xlist, Delta3, label="L={}".format(L))
# plt.plot(ETlist, Delta2, label="L={}".format(L))
plt.xlim(0, max(xlist))
plt.xlabel(r"$1/E_T^3$")
plt.ylabel(r"$\Lambda$")
plt.savefig("E0pert3.pdf")
plt.clf()


xlist = 1/ETlist**2
plt.plot(xlist, Delta4, label="L={}".format(L))
# plt.plot(ETlist, Delta2, label="L={}".format(L))
plt.xlim(0, max(xlist))
plt.title("L={}".format(L))
plt.xlabel(r"$1/E_T^2$")
plt.ylabel(r"$\Lambda$")
plt.savefig("E0pert4_L={}.pdf".format(L))
plt.clf()
