import scipy
import json
import dataset
import datetime

class Database():
    def __init__(self, dbname="../spectra.db", tablename="spectra"):
        self.db = dataset.connect('sqlite:///'+dbname)
        self.table=self.db[tablename]

    def insert(self, k, L, Emax, ren, g, spec, basisSize, neigs, m, dual, finiteVolCouplings, Er=None, cutoff=None, version=None):
        L = float(L)
        k = int(k)
        Emax = float(Emax)
        if ren!="raw" and ren!="renlocal" and ren!="rensubl":
            raise Exception()
        g = float(g)
        m = float(m)
        basisSize=int(basisSize)
        neigs=int(neigs)

        self.table.insert(dict(date=datetime.datetime.now(), k=k, L=L, m=m, Emax=Emax, ren=ren, g=g, finiteVolCouplings=finiteVolCouplings, \
                                dual=dual, spec=json.dumps(spec.tolist()), basisSize=basisSize, Er=Er, cutoff=cutoff, version=version, neigs=neigs))
