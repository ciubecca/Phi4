import numpy as np
from numpy import array
from counterterms import *

L = 5
ETmax = 20

ETlist = np.linspace(10, ETmax, 10)

c = exactct(L, ETmax)

y = array([c.ct3(ET) for ET in ETlist])
plt.plot
