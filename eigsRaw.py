import statefuncs
from  profile_support import *
import phi4
import renorm
import sys
import scipy
import numpy as np
import math
import database
from numpy import pi
from matplotlib import rc
import matplotlib.pyplot as plt

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
params = {'legend.fontsize': 8, 'lines.markersize':2.5, 'lines.marker':"o"}
plt.rcParams.update(params)
plt.style.use('ggplot')

color = {1:["b","g"], -1:["r","c"]}

m = 1
# Number of eigenvalues to compute per sector
neigs = 2


argv = sys.argv
if len(argv) < 4:
    print(argv[0], "<L> <ET> <g>")
    sys.exit(-1)

L = float(argv[1])
ET = float(argv[2])
g = float(argv[3])

theta = [pi*x for x in np.linspace(0.,0.8,20.)]
glist = [g*np.exp(1J*x) for x in theta]

print(glist)

print("L, ET", L, ET)

def main():

    eigs = {}

    for k in (-1,1):
        a = phi4.Phi4(m, L, k)
        a.buildBasis(Emax=ET)

        # Compute the potential matrices in the low-energy space below ET
        a.computePotential()

        print("Full basis size: ", a.basis.size)

        a.setglist(glist=glist)

        a.computeEigval(ET, "raw", neigs=neigs)
        eigs[k] = np.array([a.eigenvalues[g]["raw"] for g in glist]).transpose()
        # print("Raw vacuum:", eigs[k][0])

    spec = {}
    spec[1] = eigs[1][1:,:]-eigs[1][0,:]
    spec[-1] = eigs[-1]-eigs[1][0,:]

    for k in (-1,1):
        for i in range(neigs):
            if i==0:
                label1='Re, k={}'.format(k)
                label2='Im, k={}'.format(k)
            else:
                label1=None
                label2=None
            plt.plot(theta, eigs[k][i].real, label=label1, color='r')
            plt.plot(theta, eigs[k][i].imag, label=label2, color='b')

        plt.xlabel(r"$\theta$")
        plt.legend()
        plt.savefig("sp_{}.pdf".format(k))
        plt.clf()

    # for k in (-1,1):
        # for i in range(spec[k].shape[0]):
            # if i==0:
                # label1='k={}, Re'.format(k)
                # label2='k={}, Im'.format(k)
            # else:
                # label1=None
                # label2=None
            # plt.plot(theta, spec[k][i].real, label=label1, color=color[k][0])
            # plt.plot(theta, spec[k][i].imag, label=label2, color=color[k][1])

    k = -1
    label1='k={}, Re'.format(k)
    label2='k={}, Im'.format(k)
    plt.plot(theta, spec[k][0].real, label=label1, color=color[k][0])
    plt.plot(theta, spec[k][0].imag, label=label2, color=color[k][1])
    plt.xlabel(r"$\theta$")
    plt.ylabel(r"$m_{ph}$")
    plt.title(r"$\vert g\vert$={}, L={}, $E_T$={}".format(g, L, ET))
    plt.legend()
    plt.savefig("mass.pdf")
    plt.clf()


main()
