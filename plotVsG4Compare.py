import sys
import matplotlib.pyplot as plt
import scipy
import scipy.optimize
import math
import json
from scipy import pi, log
from matplotlib import rc
import database
import finiteVolH

version = "v3_5-LV"

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

def main(argv):

    if len(argv) < 5:
        print argv[0], "<L> <Emax> <mindim> <maxdim>"
        return -1

    L = float(argv[1])
    Emax = float(argv[2])
    mindim = int(argv[3])
    maxdim = int(argv[4])

    # Hardcoded parameters
    neigs = 5

    plt.figure(1) # E_0
    plt.clf()

    plt.figure(2) # E_i-E_0
    plt.clf()

    params = {'legend.fontsize': 8}
    plt.rcParams.update(params)

    xDirect = scipy.linspace(0., 3., num=16, endpoint=True)
    print xDirect

    xDual = [2.262,2.275]+list(scipy.linspace(2.3, 3.0, 15))
    print xDual

    db = database.Database('spectra.db')

####################   DUAL BASIS   #######################
    oddColor = 'red'
    evenColor = 'blue'

    vacuum = { }
    evenSpectrum = { }
    oddSpectrum = { }

    t1 = { }
    for ren in ("raw","renlocal","rensubl"):
        t1[ren] = { }
        for k in (-1,1):
            t1[ren][k] = [e for e in db.table.find(ren=ren, k=k, L=L, finiteVolCouplings=True, dual=True, version=version)]

    t2 = { }
    for ren in ("raw","renlocal","rensubl"):
        t2[ren] = { }
        for k in (-1,1):
            t2[ren][k] = []
            for g in xDual:
                listRes = []
                # Find all entries with given g4
                for e in t1[ren][k]:
                    if abs(e['g']-g)<10.**(-12.) and mindim<e['basisSize']<maxdim:
                        listRes.append((e['basisSize'], json.loads(e['spec'])))
                if listRes == []:
                    print 'Dual spectrum not found: g=', g, ', k=', k
                    return -1
                # Select the entry with highest basisSize
                e = max(listRes, key = lambda p:p[0])
                t2[ren][k].append(e[1:])
    
    xList = scipy.linspace(2.262, 3, 100)
    a = finiteVolH.FiniteVolH(L, 1)
    muList = []
    E0List = []
    for g in xList:
        g0, g2, g4 = a.dualCouplings(g)
        mu = a.mu
        muList.append(mu)
        E0List.append(L*(g0-(mu**2./2.+g2)**2./(4.*g4)))
    E0List = scipy.array(E0List)
    muList = scipy.array(muList)

    for ren in ("raw","renlocal","rensubl"):
        vacuum[ren] = scipy.array([e[0][0] for e in t2[ren][1]])
        evenSpectrum[ren] = [ scipy.array([e[0][i] for e in t2[ren][1]]) - vacuum[ren] for i in range(1,neigs) ]
        oddSpectrum[ren] = [ scipy.array([e[0][i] for e in t2[ren][-1]]) - vacuum[ren] for i in range(neigs) ]

    def label(ren, maxdim, k=None):
        p = ""
        if k==1:
            p = r", $Z_2=+$"
        elif k==-1:
            p = r"$, Z_2=-$"
        return ren+", basisSize="+str(maxdim)+p


    plt.figure(1)
    data = vacuum["rensubl"]
    plt.plot(xDual, data, linewidth=0.5, linestyle='-', marker='x', markersize=2, color=evenColor,markeredgecolor=evenColor)

    plt.plot(xList, E0List, 'k--')

    plt.figure(2)
    for i in range(neigs):
        l = None
        data = oddSpectrum["rensubl"][i]
        if i==0:
            l = r"$Z_2 = -$"
        plt.plot(xDual, data, linewidth=0.5, linestyle='-', marker='x', markersize=2, color=oddColor,markeredgecolor=oddColor, label=l)


    for i in range(neigs-1):
        l = None
        data = evenSpectrum["rensubl"][i]
        if i==0:
            l = r"$Z_2 = +$"
        plt.plot(xDual, data, linewidth=0.5, linestyle='-', marker='x', markersize=2, color=evenColor,markeredgecolor=evenColor, label=l)

    plt.plot(xList, muList, 'k--')

    # PLOT ERROR BANDS
    plt.figure(1)
    upper = vacuum["rensubl"]
    lower = vacuum["renlocal"]
    plt.fill_between(xDual, lower, upper, facecolor=evenColor, alpha=0.2, edgecolor=evenColor, linewidth=0.0)

    plt.figure(2)
    for i in range(neigs):
        upper = oddSpectrum["rensubl"][i]
        lower = oddSpectrum["renlocal"][i]
        plt.fill_between(xDual, lower, upper, facecolor=oddColor, alpha=0.2, edgecolor=oddColor, linewidth=0.0)

    for i in range(neigs-1):
        upper = evenSpectrum["rensubl"][i]
        lower = evenSpectrum["renlocal"][i]
        plt.fill_between(xDual, lower, upper, facecolor=evenColor, alpha=0.2, edgecolor=evenColor, linewidth=0.0)

########  NORMAL BASIS  ###########
    oddColor = 'green'
    evenColor = 'purple'

    vacuum = { }
    evenSpectrum = { }
    oddSpectrum = { }

    t1 = { }
    for ren in ("raw","renlocal","rensubl"):
        t1[ren] = { }
        for k in (-1,1):
            t1[ren][k] = [e for e in db.table.find(Emax=Emax, ren=ren, k=k, L=L, finiteVolCouplings=True, dual=False, version=version)]

    t2 = { }
    for ren in ("raw","renlocal","rensubl"):
        t2[ren] = { }
        for k in (-1,1):
            t2[ren][k] = []
            for g in xDirect:
                found = False
                for e in t1[ren][k]:
                    if abs(e['g']-g) < 10.**(-12.):
                        t2[ren][k].append(json.loads(e['spec']))
                        found = True
                        break
                if not found:
                    print 'Normal spectrum not found: g=', g, ', ren=', ren, ', k=', k
                    return -1

    for ren in ("raw","renlocal","rensubl"):
        vacuum[ren] = scipy.array([e[0] for e in t2[ren][1]])
        evenSpectrum[ren] = [ scipy.array([e[i] for e in t2[ren][1]]) - vacuum[ren] for i in range(1,neigs) ]
        oddSpectrum[ren] = [ scipy.array([e[i] for e in t2[ren][-1]]) - vacuum[ren] for i in range(neigs) ]


    plt.figure(1)
    data = vacuum["rensubl"]
    plt.plot(xDirect, data, linewidth=.5, dashes = [3,2], marker='.', linestyle='--', markersize=2, color=evenColor, markeredgecolor=evenColor)

    plt.figure(2)
    for i in range(neigs):
        l = None
        data = oddSpectrum["rensubl"][i]
        if i==0:
            l=r"$Z_2 = -$"
        plt.plot(xDirect, data, linewidth=.5, dashes = [3,2], marker='.', linestyle='--', markersize=2, color=oddColor, markeredgecolor=oddColor, label=l)

    for i in range(neigs-1):
        l = None
        data = evenSpectrum["rensubl"][i]
        if i==0:
            l=r"$Z_2 = +$"
        plt.plot(xDirect, data, linewidth=.5, dashes = [3,2], marker='.', linestyle='--', markersize=2, color=evenColor, markeredgecolor=evenColor, label=l)


    # PLOT ERROR BANDS
    plt.figure(1)
    upper = vacuum["rensubl"]
    lower = vacuum["renlocal"]
    plt.fill_between(xDirect, lower, upper, facecolor=evenColor, alpha=0.2, edgecolor=evenColor, linewidth=0.0)

    plt.figure(2)
    for i in range(neigs):
        upper = oddSpectrum["rensubl"][i]
        lower = oddSpectrum["renlocal"][i]
        plt.fill_between(xDirect, lower, upper, facecolor=oddColor, alpha=0.2, edgecolor=oddColor, linewidth=0.0)
    for i in range(neigs-1):
        upper = evenSpectrum["rensubl"][i]
        lower = evenSpectrum["renlocal"][i]
        plt.fill_between(xDirect, lower, upper, facecolor=evenColor, alpha=0.2, edgecolor=evenColor, linewidth=0.0)

#############################################################################


    plt.figure(1)
    plt.gcf().set_size_inches(4,3)
    plt.subplots_adjust(bottom=0.16,left=0.18,right=0.95)
    plt.ylim(-2.,0.01)
    plt.xlim(min(xDirect)-0.01, max(xDirect)+0.01)
    plt.title(r"$m=1, \quad L=$"+str(L), fontsize=10)
    plt.xlabel(r"$g$")
    plt.ylabel(r"$\mathcal{E}_0$")
    plt.legend(loc='lower left', prop={'size':6})
    plt.savefig("fig_vacCompare_L="+str(L)+".pdf", prop={'size':10})

    plt.figure(2)
    plt.gcf().set_size_inches(4,3)
    plt.subplots_adjust(bottom=0.16,left=0.15,right=0.95)
    plt.ylim(0.,4.8)
    plt.xlim(min(xDirect)-0.01, max(xDirect)+0.01)
    plt.title(r"$m=1, \quad L=$"+str(L), fontsize=10)
    plt.xlabel(r"$g$")
    plt.ylabel(r"$\mathcal{E}_I-\mathcal{E}_0$")
    plt.legend(loc='lower left', prop={'size':6})
    plt.savefig("fig_specCompare_L="+str(L)+".pdf", prop={'size':10})


if __name__ == "__main__":
    main(sys.argv)
