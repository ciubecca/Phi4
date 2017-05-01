import sys
import scipy
from scipy.optimize import curve_fit
import math
from scipy import pi, log, log10, array, sqrt, stats, exp, e
import database
from sys import exit
import numpy as np
from numpy import testing
from sklearn.linear_model import Ridge, LinearRegression
from scipy.special import kn
from itertools import combinations

LList = [5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10]

# Minimum g for which we use nonlinear fit
glim = 1.

# Fitting vs ET, number of lowest values of ET to potentially exclude
nmax = 5
# Maximum number of points to exclude at low ET
pmax = 2


ETmax = {}
ETmin = {}

ETmax["rentails"] = {5:31, 5.5:29, 6:27.5, 6.5:26.5, 7:25, 7.5:24, 8:23,
        8.5:22, 9:21, 9.5:20.5, 10:20}

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

# TODO Probably need to do the fit on the subtracted spectrum directly
class Extrapolator():
    """ Extrapolate vs ET """

    def __init__(self, db, k, L, g, ren="rentails"):
        ETMin = ETmin[ren][L]
        ETMax = ETmax[ren][L]

        self.g = g
        self.L = L
        self.k = k

        mult = 1/step[ren]

        if ren != "rentails":
            self.ETlist = np.array([ETMin + n for n in range(ETMax-ETMin+1) if
                ETMin+n not in missing[L]])
        else:
            self.ETlist = scipy.linspace(ETMin, ETMax, (ETMax-ETMin)*2+1)

        self.spectrum = np.array([np.array(db.getEigs(k, ren, g, L, ET))
            for ET in self.ETlist])

    def train(self, neigs=1):

        self.neigs = neigs

        # Use just linear fit because of fluctuations
        if self.g < glim:
            self.featureVec = lambda ET: [1/ET**3]
        else:
            self.featureVec = lambda ET: [1/ET**3, 1/ET**4]


        # Position of the highest ET to exclude
        self.models = []

        # Number of points to exclude
        for n in range(neigs):
            self.models.append([])

            for m in range(0,pmax+1):
                for nlist in combinations(range(nmax+1), m):
                    mask = np.ones(len(self.ETlist), dtype=bool)
                    mask[list(nlist)] = False
                    data = self.spectrum[mask, n]
                    xlist = self.ETlist[mask]
                    X = np.array(self.featureVec(xlist)).transpose()
                    self.models[n].append(LinearRegression().fit(X, data))

    def predict(self, x):
        x = np.array(x)
        N = len(self.models[0])
        return np.array([sum(self.models[m][n].predict(np.array(self.featureVec(x)).\
                        transpose()) for n in range(N))/N for m in range(self.neigs)])

    def asymValue(self):
        N = len(self.models[0])
        ints = np.array([[self.models[m][n].intercept_ for n in range(N)]
            for m in range(self.neigs)])
        return np.mean(ints, axis=1)

    def asymErr(self):
        # Number of models
        N = len(self.models[0])


        # Intercepts
        ints = np.array([[self.models[m][n].intercept_ for n in range(N)]
            for m in range(self.neigs)])

        asymVal = self.asymValue()

        # Residuals
        predictions = np.array([[self.models[m][n].predict(
            np.array(self.featureVec(self.ETlist)).transpose()) for n in range(N)]
            for m in range(self.neigs)])
        residuals = self.spectrum[:, np.newaxis, :self.neigs].transpose() - predictions

        # Intercepts plus max residuals of the last 10 points
        def positive(x):
            return (x>0)*x
        upperBounds = ints + np.amax(positive(residuals)[:,:,-10:], axis=2)
        # Select the highest upper bound among all the models
        upperBound = np.amax(upperBounds, axis=1)

        # Intercepts plus min residuals of the last 10 points
        def negative(x):
            return (x<0)*x
        lowerBounds = ints + np.amin(negative(residuals)[:,:,-10:], axis=2)
        # Select the highest upper bound among all the models
        lowerBound = np.amin(lowerBounds, axis=1)

        testing.assert_array_less(lowerBound, asymVal)
        testing.assert_array_less(asymVal, upperBound)

        ret = np.array([asymVal-lowerBound,upperBound-asymVal]).transpose()

        # Return two-dim array (lowerErr, upperErr)
        return ret

def Massfun(L, m, b, c):
    return m*(1 + b*kn(1, m*L) + c/(L*m)**(3/2)*e**(-m*L))
fmassStr = r"$m_{ph}(1 + b  K_1(m_{ph} L) + c/(L m_{ph})^{3/2} e^{-m_{ph}L} )$"

def Massfun2(L, m, b):
    return m*(1 + b*kn(1, m*L))
fmassStr2 = r"$m_{ph}(1 + b  K_1(m_{ph} L))$"

def Lambdafun(L, a, m):
    return a - m/(pi*L)*kn(1, m*L)
fvacStr = r"$\Lambda - \frac{m_{ph}}{\pi L} K_1(m_{ph} L)$"

method = 'dogbox'

# XXX maybe we should fit the subtracted spectrum directly!
class ExtrvsL():
    def __init__(self, db, g):

        self.g = g

        self.LambdaInf = np.zeros(len(LList))
        self.LambdaErr = np.zeros((2, len(LList)))
        self.MassInf = np.zeros(len(LList))
        self.MassErr = np.zeros((2, len(LList)))

        for i,L in enumerate(LList):
            e = {}
            e[1] = Extrapolator(db, 1, L, g)
            e[-1] = Extrapolator(db, -1, L, g)
            e[1].train()
            e[-1].train()

            self.LambdaInf[i] = e[1].asymValue()[0]/L
            self.LambdaErr[:,i] = e[1].asymErr()[0]/L

            self.MassInf[i] = e[-1].asymValue()[0]-e[1].asymValue()[0]
            # XXX Check
            self.MassErr[:,i] = np.amax([e[-1].asymErr()[0],
                e[1].asymErr()[0][::-1]], axis=0)
            # self.MassErr[:,i] = e[-1].asymErr()+e[1].asymErr()[::-1]

    def train(self, nparam=3):
        """ if nparam=2, use only 2 free coefficients for fitting the mass.
        This should be better in the critical region, where S=-1 """

        self.popt = {k: [None, None] for k in (-1,1)}
        self.msg = {}
        self.coefs = {k: [] for k in (-1,1)}
        self.errs = {k: [] for k in (-1,1)}

        popt = self.popt
        coefs = self.coefs
        errs = self.errs

        # Upper or lower values
        try:
            for n in (0,1):
                y = self.MassInf -(-1)**n*self.MassErr[n]
                if nparam==2:
                    self.mfun = Massfun2
                else:
                    self.mfun = Massfun
                # Weigh points by error bars
                popt[-1][n], pcov = curve_fit(self.mfun, LList, y.ravel(),
                        method=method, sigma=self.MassErr[n])

            # Estimate of physical mass
            # mph = (+popt[-1][1][0])/2
            # self.LambdaFixedM = lambda L,a: Lambdafun(L, a, mph)

            # Estimated bounds on the physical mass
            bounds = ([-np.inf, popt[-1][0][0]], [np.inf, popt[-1][1][0]])

            for n in (0,1):
                # NOTE Fix Mph instead of fitting it
                y = self.LambdaInf -(-1)**n*self.LambdaErr[n]
                # Weigh points by error bars
                popt[1][n], pcov = curve_fit(Lambdafun, LList, y.ravel(),
                        bounds=bounds, method=method, sigma=self.LambdaErr[n])

        except RuntimeError as e:
            print("Exception for g={}".format(self.g))
            raise e

        for k in (-1,1):
            for i in range(len(popt[k][0])):
                x = popt[k][0][i], popt[k][1][i]
                coefs[k].append((x[0]+x[1])/2)
                errs[k].append(abs(x[0]-x[1])/2)

        self.msg[1] = [
            r"$\Lambda = {:.7f} \pm {:.7f}$".format(coefs[1][0], errs[1][0])
            ,r"$m_{{ph}} = {:.7f} \pm {:.7f}$".format(coefs[1][1], errs[1][1])
        ]

        self.msg[-1] = [
            r"$m_{{ph}} = {:.7f} \pm {:.7f}$".format(coefs[-1][0], errs[-1][0])
            , r"$b = {:.7f} \pm {:.7f}$".format(coefs[-1][1], errs[-1][1])
        ]
        if nparam==3:
            self.msg[-1].append(
                    r"$c = {:.7f} \pm {:.7f}$".format(coefs[-1][2], errs[-1][2])
                    )

    def predict(self, k, x):
        if k==1:
            fun = Lambdafun
        else:
            fun = self.mfun

        return (fun(x, *self.popt[k][1])+fun(x, *self.popt[k][0]))/2

    def asymValue(self):
        return {k: self.coefs[k][0] for k in (-1,1)}

    def asymErr(self):
        return {k: self.errs[k][0] for k in (-1,1)}
