import sys
import scipy
from scipy.optimize import curve_fit
import math
from scipy import pi, log, log10, array, sqrt, stats
import database
from sys import exit
import numpy as np
from sklearn.linear_model import Ridge


ETmax = {5:32, 5.5:30, 6:28, 6.5:26.5, 7:25, 7.5:23.5, 8:23,
        8.5:21.5, 9:21, 9.5:20, 10:20}
ETmax[10] = 19.5
ETmin = {5:15, 5.5:15, 6:10, 6.5:10, 7:10, 7.5:10, 8:10,
        8.5:10, 9:10, 9.5:10, 10:10}

def stdWeights(ET):
    return 1

absolute_sigma=False
errcoeff = 3.


def stdFeatureVec(ET):
    return [10/ET**2, 1/ET**3, log(ET**(1/ET**2))]

def stdFeatureVec(ET):
    return [1/ET**3, 1/ET**4]

class Extrapolator():

    def __init__(self, db, k, L, g, ren="rentails"):
        self.ETlist = scipy.linspace(ETmin[L], ETmax[L], (ETmax[L]-ETmin[L])*2+1)
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

    def predict(self, x):
        x = np.array(x)
        return self.model.predict(np.array(self.featureVec(x)).transpose())

    def asymptoticValue(self):
        return self.model.intercept_
