import numpy as np
from numpy import array
from counterterms import *
from paramplots import *

L = 5
m = 1
ETmin = 5
ETmax = 15

ETlist = np.linspace(ETmin, ETmax, 10)

c = exactct(L, ETmax)

y = array([c.ct3(ET) for ET in ETlist])
plt.plot(ETlist, y, label="exact")
y = array([ct0ET3(ET, m) for ET in ETlist])
plt.plot(ETlist, y)

plt.legend()
plt.savefig("test.pdf")
