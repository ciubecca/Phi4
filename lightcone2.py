import statefuncs
import phi4
import renorm
import sys
from scipy import sqrt, array
import scipy
import math
from matplotlib import rc
from cycler import cycler
import matplotlib.pyplot as plt
from bisect import bisect_left

output = "png"


# saveondb = False
m = 1
# Number of eigenvalues to compute per sector
neigs = 1

glist = [0., 0.0151515, 0.030303, 0.0454545, 0.0606061, 0.0757576, 0.0909091, 0.106061, 0.121212, 0.136364, 0.151515, 0.166667, 0.181818, 0.19697, 0.212121, 0.227273, 0.242424, 0.257576, 0.272727, 0.287879, 0.30303, 0.318182, 0.333333, 0.348485, 0.363636, 0.378788, 0.393939, 0.409091, 0.424242, 0.439394, 0.454545, 0.469697, 0.484848, 0.5, 0.515152, 0.530303, 0.545455, 0.560606, 0.575758, 0.590909, 0.606061, 0.621212, 0.636364, 0.651515, 0.666667, 0.681818, 0.69697, 0.712121, 0.727273, 0.742424, 0.757576, 0.772727, 0.787879, 0.80303, 0.818182, 0.833333, 0.848485, 0.863636, 0.878788, 0.893939, 0.909091, 0.924242, 0.939394, 0.954545, 0.969697, 0.984848, 1., 1.01515, 1.0303, 1.04545, 1.06061, 1.07576, 1.09091, 1.10606, 1.12121, 1.13636, 1.15152, 1.16667, 1.18182, 1.19697, 1.21212, 1.22727, 1.24242, 1.25758, 1.27273, 1.28788, 1.30303, 1.31818, 1.33333, 1.34848, 1.36364, 1.37879, 1.39394, 1.40909, 1.42424, 1.43939, 1.45455, 1.4697, 1.48485, 1.5]


msqlist = [1., 0.999666, 0.998697, 0.997136, 0.99502, 0.992378, 0.989235, 0.985612, 0.981528, 0.976999, 0.972038, 0.966657, 0.960867, 0.954678, 0.948096, 0.941131, 0.933788, 0.926073, 0.917991, 0.909546, 0.900743, 0.891584, 0.882074, 0.872215, 0.862009, 0.851459, 0.840565, 0.82933, 0.817756, 0.805842, 0.79359, 0.781, 0.768073, 0.754809, 0.741209, 0.727272, 0.712998, 0.698386, 0.683437, 0.668148, 0.652521, 0.636553, 0.620244, 0.603592, 0.586597, 0.569256, 0.551569, 0.533533, 0.515147, 0.496409, 0.477317, 0.457868, 0.438062, 0.417894, 0.397363, 0.376466, 0.3552, 0.333562, 0.31155, 0.28916, 0.266389, 0.243233, 0.219689, 0.195753, 0.171422, 0.14669, 0.121554, 0.0960092, 0.070051, 0.0436746, 0.0168749, -0.0103535, -0.0380161, -0.0661186, -0.0946673, -0.123668, -0.153129, -0.183055, -0.213455, -0.244337, -0.275708, -0.307578, -0.339955, -0.372849, -0.406271, -0.440233, -0.474745, -0.50982, -0.545474, -0.58172, -0.618576, -0.656061, -0.694194, -0.733002, -0.772527, -0.883792, -1.05411, -1.22565, -1.39842, -1.57246]


msqlistren = [1., 0.999665, 0.998694, 0.997129, 0.995005, 0.992352, 0.989192, 0.985543, 0.981424, 0.976849, 0.971816, 0.96635, 0.960451, 0.954122, 0.947364, 0.940176, 0.932556, 0.924497, 0.916103, 0.90717, 0.89777, 0.887885, 0.877774, 0.867248, 0.856304, 0.844941, 0.833156, 0.820946, 0.808308, 0.79524, 0.781737, 0.767796, 0.753413, 0.738583, 0.723303, 0.707567, 0.691371, 0.67471, 0.657578, 0.639971, 0.621882, 0.603306, 0.584237, 0.564669, 0.544595, 0.52401, 0.502906,
        0.481276, 0.459114, 0.436412, 0.413163, 0.38936, 0.364994, 0.340057, 0.314541, 0.288438, 0.261738, 0.234434, 0.206515, 0.177972, 0.148795, 0.118975, 0.0885005, 0.0573611, 0.0255458, -0.00695687, -0.0401586, -0.0740716, -0.108708, -0.144082, -0.180206, -0.217095, -0.254762, -0.293224, -0.332495, -0.372593, -0.413534, -0.455336, -0.498018, -0.5416, -0.586102, -0.631546, -0.677955, -0.725352, -0.773765, -0.82322, -0.873747, -0.925376, -0.978142, -1.03208, -1.08723,
        -1.14364, -1.20136, -1.26044, -1.32096, -1.45786, -1.65792, -1.95263, -2.27165, -2.65325]

glist = glist[:80]
msqlist = msqlist[:80]
msqlistren = msqlistren[:80]

gETlist = scipy.linspace(0,3.,31)

print("gETlist", gETlist)


argv = sys.argv
if len(argv) < 3:
    print(argv[0], "<L> <ET>")
    sys.exit(-1)

L = float(argv[1])
ET = float(argv[2])


a = phi4.Phi4(m, L)
a.buildBasis(Emax=ET)


print("Full basis size: ", a.basis[1].size)


a.computePotential(k=1)
a.computePotential(k=-1)

phi2 = a.V[1][2].M/L

mET = []
mETren = []
VEV = []
VEVren = []

for g in gETlist:

    print("g ET=", g)
    print("Computing raw eigenvalues")

    a.setCouplings(g4=g)
    a.computeEigval(k=1, ET=ET, ren="raw", neigs=neigs)
    a.computeEigval(k=-1, ET=ET, ren="raw", neigs=neigs)
    v0 = scipy.array(a.eigenvectors["raw"][1][0])

    VEV.append(phi2.dot(v0).dot(v0))
    mET.append(a.spectrum(-1,"raw")[0])

    eps = a.eigenvalues["raw"][1][0]
    a.computeEigval(k=1, ET=ET, ren="renloc", EL=ET, eps=eps, neigs=neigs)
    a.computeEigval(k=-1, ET=ET, ren="renloc", EL=ET, eps=eps, neigs=neigs)
    v0 = scipy.array(a.eigenvectors["renloc"][1][0])
    VEVren.append(phi2.dot(v0).dot(v0))
    mETren.append(a.spectrum(-1,"renloc")[0])



params = {'legend.fontsize': 8}
plt.rcParams.update(params)

plt.rc('axes', prop_cycle=(cycler('color', ['r', 'g', 'b', 'y']) +
    cycler('linestyle', ['-', '--', ':', '-.'])) +
    cycler('label', ["ET","LT","LT extrapolated", "ET ren"]))

mLC = []
mLCextr = []
mLCren = []
mLCextrren = []


for i,g in enumerate(gETlist):

    msqLC = m**2+12*g*VEV[i]
    gLC = g/msqLC

    j = bisect_left(glist, gLC)
    mLC.append(sqrt(msqlist[j]*msqLC))
    mLCextr.append(sqrt(msqlistren[j]*msqLC))


for i,g in enumerate(gETlist):

    msqLC = m**2+12*g*VEVren[i]
    gLC = g/msqLC

    j = bisect_left(glist, gLC)
    mLCren.append(sqrt(msqlist[j]*msqLC))
    mLCextrren.append(sqrt(msqlistren[j]*msqLC))


xlist = gETlist

plt.figure(1)

for ylist in (mET, mLC, mLCextr):
    plt.plot(xlist, ylist)

plt.figure(2)

for ylist in (mETren, mLCren, mLCextrren):
    plt.plot(xlist, ylist)

title = r"$L$={:.1f}, $E_T$={:.1f}".format(L,ET)
fname = "comparison2LT_L={0:.1f}_ET={1:.1f}.{2}".format(L,ET,output)
loc = "lower right"

plt.figure(1, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
plt.title(title)
plt.ylim(0,1.1)
plt.xlabel(r"$g$")
plt.ylabel(r"$m$")
plt.legend(loc=loc)

plt.savefig(fname)




title = r"$L$={:.1f}, $E_T$={:.1f} ren".format(L,ET)
fname = "comparison2LTren_L={0:.1f}_ET={1:.1f}.{2}".format(L,ET,output)
loc = "lower right"

plt.figure(2, figsize=(4., 2.5), dpi=300, facecolor='w', edgecolor='w')
plt.title(title)
plt.ylim(0,1.1)
plt.xlabel(r"$g$")
plt.ylabel(r"$m$")
plt.legend(loc=loc)

plt.savefig(fname)
