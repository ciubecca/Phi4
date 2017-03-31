import sys
import scipy
from scipy.optimize import curve_fit
import math
from scipy import pi, log, log10, array, sqrt, stats
import database
from sys import exit
import numpy as np
from sklearn.linear_model import Ridge

ETmax = {}
ETmin = {}

ETmax["rentails"] = {5:32, 5.5:30, 6:28, 6.5:26.5, 7:25, 7.5:24, 8:23,
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
        xlist = self.ETlist
# Divide by L
        data = self.spectrum

        if weights==None:
            weights = stdWeights
        if featureVec==None:
            self.featureVec = stdFeatureVec
        else:
            self.featureVec = featureVec

        X = np.array(self.featureVec(self.ETlist)).transpose()
        self.model = Ridge(alpha=alpha, normalize=True)
        self.model.fit(X, data, sample_weight=1/weights(xlist))
        print("L={}, k={} fit:".format(self.L,self.k)+str(self.model.coef_))

    def predict(self, x):
        x = np.array(x)
        return self.model.predict(np.array(self.featureVec(x)).transpose())

    def asymValue(self):
        return self.model.intercept_

    def asymErr(self):
        return max(abs(self.predict(self.ETlist)-self.spectrum)[-10:])
