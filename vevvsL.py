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
from sklearn.linear_model import LinearRegression

# Inverse power of L dependence of the VEV
power = 1

gLClist = np.array([0., 0.0120368, 0.0240735, 0.0361103, 0.048147, 0.0601838, 0.0722205, 0.0842573, 0.096294, 0.108331, 0.120368, 0.132404, 0.144441, 0.156478, 0.168515, 0.180551, 0.192588, 0.204625, 0.216662, 0.228698, 0.240735, 0.252772, 0.264809, 0.276845, 0.288882, 0.300919, 0.312956, 0.324992, 0.337029, 0.349066, 0.361103, 0.373139, 0.385176, 0.397213, 0.40925, 0.421286, 0.433323, 0.44536, 0.457397, 0.469433, 0.48147, 0.493507, 0.505544, 0.51758, 0.529617, 0.541654, 0.553691, 0.565727, 0.577764,
0.589801, 0.601838, 0.613874, 0.625911, 0.637948, 0.649985, 0.662021, 0.674058, 0.686095, 0.698132, 0.710168, 0.722205, 0.734242, 0.746279, 0.758315, 0.770352, 0.782389, 0.794426, 0.806462, 0.818499, 0.830536, 0.842573, 0.854609, 0.866646, 0.878683, 0.89072, 0.902757, 0.914793, 0.92683, 0.938867, 0.950904, 0.96294, 0.974977, 0.987014, 0.999051, 1.01109, 1.02312, 1.03516, 1.0472, 1.05923, 1.07127, 1.08331, 1.09534, 1.10738, 1.11942, 1.13145, 1.14349, 1.15553, 1.16757, 1.1796, 1.19164, 1.20368,
1.21571, 1.22775, 1.23979, 1.25182, 1.26386, 1.2759, 1.28793, 1.29997, 1.31201, 1.32404, 1.33608, 1.34812, 1.36015, 1.37219, 1.38423, 1.39626, 1.4083, 1.42034, 1.43237, 1.44441, 1.45645, 1.46848, 1.48052, 1.49256, 1.50459, 1.51663, 1.52867, 1.5407, 1.55274, 1.56478, 1.57681, 1.58885, 1.60089, 1.61292, 1.62496, 1.637, 1.64904, 1.66107, 1.67311, 1.68515, 1.69718, 1.70922, 1.72126, 1.73329, 1.74533, 1.75737, 1.7694, 1.78144, 1.79348, 1.80551, 1.81755, 1.82959, 1.84162, 1.85366, 1.8657,
1.87773, 1.88977, 1.90181, 1.91384, 1.92588, 1.93792, 1.94995, 1.96199, 1.97403, 1.98606, 1.9981, 2.01014, 2.02217, 2.03421, 2.04625, 2.05828, 2.07032, 2.08236, 2.0944])

msqLClist2 = np.array([1., 0.999787, 0.999167, 0.998162, 0.996792, 0.995073, 0.993018, 0.990642, 0.987954, 0.984962, 0.981676, 0.9781, 0.974241, 0.970104, 0.965701, 0.96102, 0.95607, 0.950869, 0.945388, 0.939641, 0.93366, 0.927384, 0.920888, 0.914078, 0.906991, 0.899713, 0.892074, 0.88427, 0.876059, 0.867718, 0.858909, 0.850018, 0.840578, 0.831118, 0.821002, 0.81095, 0.800617, 0.789999, 0.779093, 0.767896, 0.756405, 0.744616, 0.732526, 0.720131, 0.707427, 0.694409, 0.681075, 0.667419, 0.653437,
0.639124, 0.624475, 0.609486, 0.594152, 0.578466, 0.562424, 0.54602, 0.529248, 0.512103, 0.494577, 0.476665, 0.458361, 0.439657, 0.420546, 0.401022, 0.381077, 0.360704, 0.339894, 0.31864, 0.296933, 0.274765, 0.252127, 0.229009, 0.205404, 0.1813, 0.156687, 0.131556, 0.105896, 0.0796944, 0.0529413, 0.0256242, -0.00226926, -0.0307518, -0.0598369, -0.0895382, -0.11987, -0.150848, -0.182486, -0.214803, -0.247814, -0.281536, -0.31599, -0.351194, -0.387169, -0.423936, -0.461518,
-0.499939, -0.539223, -0.579399, -0.620494, -0.66254, -0.705567, -0.749612, -0.794712, -0.840908, -0.888246, -0.936774, -0.991022, -1.11862, -1.31174, -1.5683, -1.86286, -2.20153, -2.56296, -2.95981, -3.37537, -3.80286, -4.2454, -4.69567, -5.15299, -5.61393, -6.07406, -6.53904, -6.99382, -7.45283, -7.91143, -8.35523, -8.80465, -9.22485, -9.65204, -10.0876, -10.5113, -10.915, -11.3279, -11.7499, -12.1394, -12.5182, -12.9069, -13.3058, -13.7151, -14.1347, -14.5649, -15.0057, -15.4569,
-15.9185, -16.3906, -16.8728, -17.3652, -17.8676, -18.3799, -18.9019, -19.4336, -19.9749, -20.5258, -21.0861, -21.6559, -22.2351, -22.8236, -23.4215, -24.0285, -24.6447, -25.2699, -25.9041, -26.547, -27.1985, -27.8585, -28.5267, -29.203, -29.8872, -30.5791, -31.2785, -31.9852, -32.6991, -33.4199, -34.1475, -34.8817])


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
neigs = 2

argv = sys.argv
if len(argv) < 2:
    print(argv[0], " <g>")
    sys.exit(-1)

g = float(argv[1])

Llist = np.linspace(5,12,8)
print("Llist: ", Llist)

glist = np.array([g])

ETlist = np.array([15,18,20])
# ETlist = np.array([18])

vevren = {ET:[] for ET in ETlist}
vevraw = {ET:[] for ET in ETlist}
vacren = {ET:[] for ET in ETlist}
vacraw = {ET:[] for ET in ETlist}

for L in Llist:
    for ET in ETlist:

        a = phi4.Phi4(m, L, 1)
        a.buildBasis(Emax=ET)
        a.computePotential()
        a.setglist(glist=glist)

        a.computeEigval(ET, "raw", neigs=neigs)
        E0raw = {g: a.eigenvalues[g]["raw"][0] for g in glist}
        vacraw[ET].append(E0raw[glist[0]])
        vevraw[ET].append(a.vev[glist[0]]["raw"])

        a.computeEigval(ET, "renloc", neigs=neigs, eps=E0raw)
        E0ren = {g: a.eigenvalues[g]["renloc"][0] for g in glist}
        vacren[ET].append(E0ren[glist[0]])
        vevren[ET].append(a.vev[glist[0]]["renloc"])


f = LinearRegression()
X = (1/Llist**power).reshape(-1, 1)
y = vevren[max(ETlist)]
f.fit(X, y)

plt.figure(1)

for ET in ETlist:
    # plt.plot(Llist, vevraw[g], label="raw, g={}".format(g))
    plt.plot(1/Llist**power, vevren[ET], label="Emax={}".format(ET))
# plt.xlim(xmin, xmax)
# plt.ylim(0,1)

X = np.linspace(0, max(1/Llist**power), 100).reshape(-1, 1)
plt.plot(X, f.predict(X), label = "fit", marker=None)


# plt.figure(2)
# for g in glist:
    # plt.plot(Llist, vacraw[g], label="raw, g={}".format(g))
    # plt.plot(Llist, vacren[g], label="ren, g={}".format(g))
# plt.xlim(xmin, xmax)
# plt.ylim(0,1)



plt.figure(1)
# plt.ylim(0.8,1)
plt.xlim(0, max(1/Llist**power))
plt.xlabel(r"$1/L^{}$".format(power))
plt.ylabel(r"$\langle \phi^2 \rangle$")
plt.title(r"$g$={}".format(g))
plt.legend()
fname = ".{0}".format(output)
s = "VEVvsL_g={}".format(g)
# plt.savefig(s+fname, bbox_inches='tight')
plt.savefig(s+fname)


# plt.figure(2)
# # plt.ylim(0.8,1)
# plt.xlabel(r"$L$")
# plt.ylabel(r"$E_0$")
# plt.title(r"$g$={}, $E_T$= {}".format(g,ET))
# plt.legend()
# fname = ".{0}".format(output)
# s = "E0vsL_g={}_Emax={}".format(g,ET)
# # plt.savefig(s+fname, bbox_inches='tight')
# plt.savefig(s+fname)
