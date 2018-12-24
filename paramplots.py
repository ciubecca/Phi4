import matplotlib.pyplot as plt
from matplotlib import rc
from cycler import cycler
import matplotlib as mpl

paramd = {'lines.linestyle': ['-', '--', ':', '-.'],
        'lines.marker': ['o','v','.','-'],
        'lines.markersize': [1,2,1,1]}

# print(mpl.rcParams.keys())
plt.style.use('ggplot')

# rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)


color = {1:"b", -1:"r"}


# plt.rc('axes', prop_cycle=(cycler('color', ['r', 'g', 'b', 'y']) +
    # cycler('linestyle', ['-', '--', ':', '-.'])))


params = {'legend.fontsize': 8, 'figure.autolayout': True}
plt.rcParams.update(params)


def setparams(idx):
    for k,v in paramd.items():
        mpl.rcParams[k] = v[idx]
