import sys
import scipy
from scipy.optimize import curve_fit
import math
from scipy import pi, log, log10, array, sqrt, stats, exp
import database
from sys import exit
import numpy as np
from sklearn.linear_model import Ridge, LinearRegression
from scipy.special import kn
from itertools import combinations

LList = [5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10]
# LList = [5,6,7,8,9,10]

ETmax = {}
ETmin = {}

ETmax["rentails"] = {5:31, 5.5:29, 6:27.5, 6.5:26.5, 7:25, 7.5:24, 8:23,
        8.5:22, 9:21, 9.5:20.5, 10:20}

# ETmax["rentails"] = {5:31, 5.5:29, 6:27, 6.5:25.5, 7:24, 7.5:22, 8:22,
         # 8.5:20.5, 9:20.5, 9.5:19.5, 10:19}

ETmin["rentails"] = {5:10, 5.5:10, 6:10, 6.5:10, 7:10, 7.5:10, 8:10,
        8.5:10, 9:10, 9.5:10, 10:10}

ETmax["raw"] = {6:50, 8:38, 10:34}
ETmax["renloc"] = ETmax["raw"]

ETmin["raw"] = ETmin["rentails"]
ETmin["renloc"] = ETmin["rentails"]

step = {}
step["raw"] = 1
step["renloc"] = step["raw"]
step["rentails"] = 0.5

missing = {6:[], 8:[32], 10:[]}

def stdWeights(ET):
    return 1

absolute_sigma=False
errcoeff = 3.


def stdFeatureVec(ET):
    return [1/ET**3, 1/ET**4]

class Extrapolator():

    def __init__(self, db, k, L, g, ren="rentails"):
        ETMin = ETmin[ren][L]
        ETMax = ETmax[ren][L]

        self.L = L
        self.k = k

        mult = 1/step[ren]

        if ren != "rentails":
            self.ETlist = np.array([ETMin + n for n in range(ETMax-ETMin+1) if
                ETMin+n not in missing[L]])
        else:
            self.ETlist = scipy.linspace(ETMin, ETMax, (ETMax-ETMin)*2+1)

        self.spectrum = np.array([db.getEigs(k, ren, g, L, ET)[0]
            for ET in self.ETlist])

    def train(self, alpha, weights=None, featureVec=None):

        if weights==None:
            weights = stdWeights
        if featureVec==None:
            self.featureVec = stdFeatureVec
        else:
            self.featureVec = featureVec


        # Position of the highest ET to exclude
        nmax = 5
        self.models = []

        # Number of points to exclude
        for m in (0,1,2):
            for nlist in combinations(range(nmax+1), m):
                # print(nlist)
                mask = np.ones(len(self.ETlist), dtype=bool)
                mask[list(nlist)] = False
                data = self.spectrum[mask]
                xlist = self.ETlist[mask]
                X = np.array(self.featureVec(xlist)).transpose()
                self.models.append(Ridge(alpha=alpha, normalize=True)\
                        .fit(X, data, sample_weight=1/weights(xlist)))


    def predict(self, x):
        x = np.array(x)
        N = len(self.models)
        return sum(self.models[n].predict(np.array(self.featureVec(x)).transpose())\
                for n in range(N))/N

    def asymValue(self):
        N = len(self.models)
        ints = np.array([self.models[n].intercept_ for n in range(N)])
        return np.mean(ints)

    def asymErr(self):
        N = len(self.models)
        ints = np.array([self.models[n].intercept_ for n in range(N)])
        asymVal = self.asymValue()
        return np.array([
            max(max(asymVal-ints),
            max((self.predict(self.ETlist)-self.spectrum)[-10:])),
            max(max(ints-asymVal),
            max((self.spectrum-self.predict(self.ETlist))[-10:]))
            ])

def Massfun(L, m, b, c):
    return m + b/L*kn(1, m*L) + c*exp(-m*L)/(L)**(5/2)
fmassStr = r"$m_{ph} + \frac{b}{L} K_1(m_{ph} L) + c\,e^{-m_{ph}L}\frac{1}{L^{5/2}}$"

def Lambdafun(L, a, m):
    return a - m/(pi*L)*kn(1, m*L)
fvacStr = r"$\Lambda - \frac{m_{ph}}{\pi L} K_1(m_{ph} L)$"


class ExtrvsL():
    def __init__(self, db, g, alpha=0):

        self.LambdaInf = np.zeros(len(LList))
        self.LambdaErr = np.zeros((2, len(LList)))
        self.MassInf = np.zeros(len(LList))
        self.MassErr = np.zeros((2, len(LList)))

        for i,L in enumerate(LList):
            e = {}
            e[1] = Extrapolator(db, 1, L, g)
            e[-1] = Extrapolator(db, -1, L, g)
            e[1].train(alpha)
            e[-1].train(alpha)
            self.LambdaInf[i] = e[1].asymValue()/L
            self.LambdaErr[:,i] = e[1].asymErr()/L
            self.MassInf[i] = e[-1].asymValue()-e[1].asymValue()
            # XXX Check
            self.MassErr[:,i] = np.amax([e[-1].asymErr(),e[1].asymErr()[::-1]], axis=0)
            # MassErr[:,i] = e[-1].asymErr()+e[1].asymErr()[::-1]

    def train(self):

        self.popt = {}
        self.pcov = {}
        self.msg = {}

        popt = self.popt
        pcov = self.pcov

        popt[1], pcov[1] = curve_fit(Lambdafun, LList, self.LambdaInf.ravel())
        self.msg[1] = [
            r"$\Lambda = {:.7f} \pm {:.7f}$".format(popt[1][0],
                np.sqrt(pcov[1][0,0])),
            r"$m_{{ph}} = {:.7f} \pm {:.7f}$".format(popt[1][1],
                np.sqrt(pcov[1][1,1]))
            # ,r"$b = {:.7f} \pm {:.7f}$".format(popt[2],np.sqrt(pcov[2,2]))
        ]

        popt[-1], pcov[-1] = curve_fit(Massfun, LList, self.MassInf.ravel())
        self.msg[-1] = [
            r"$m_{{ph}} = {:.7f} \pm {:.7f}$"\
                .format(popt[-1][0], np.sqrt(pcov[-1][0,0]))
            , r"$b = {:.7f} \pm {:.7f}$".format(popt[-1][1], np.sqrt(pcov[-1][1,1]))
            , r"$c = {:.7f} \pm {:.7f}$".format(popt[-1][2], np.sqrt(pcov[-1][2,2]))
        ]

    def asymValue(self):
        return {k: self.popt[k][0] for k in (-1,1)}

    def asymErr(self):
### XXX CHECK
        return {1: np.mean(self.LambdaErr, axis=1), \
                -1: np.mean(self.MassErr, axis=1)}
