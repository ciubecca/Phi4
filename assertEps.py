import numpy as np
import sys
import scipy
import math
from scipy import pi, log, array, sqrt
import database

def main(argv):

    db = database.Database()

    for x in db.table:
        ren = x["ren"]
        if ren=="raw":
            continue
        elif ren=="renloc":
            renEps = "raw"
        elif ren=="rentails":
            renEps = "renloc"

        E = sorted(db.getEigs(x["k"],renEps,x["g"],x["L"],x["ET"]))[0]

        np.testing.assert_almost_equal(x["eps"], E)


if __name__ == "__main__":
    main(sys.argv)
