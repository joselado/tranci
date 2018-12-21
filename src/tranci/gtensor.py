# routines to compute the gtensor
import scipy.linalg as lg
import numpy as np

def get_gtensor(atom,h0):
    """Compute the gtensor according to the definition g = d2E/dBdS"""
    raise # this is not ok
    def gete(db,ds): # get Hamiltonian for infinitesimal db and ds
        out = h0 # first term
        out = out + db[0]*(atom.lx + 2.*ds[0]*atom.sx) # first component
        out = out + db[1]*(atom.ly + 2.*ds[1]*atom.sy) # second component
        out = out + db[2]*(atom.lz + 2.*ds[2]*atom.sz) # third component
        return lg.eigvalsh(out.todense())[0] # return lowest energy
#        return out # return Hamiltonian
    h = 1e-3 # infinitesimal
    out = np.zeros((3,3))  # initialize
    for i in range(3): # first loop
        for j in range(3): # second loop
            db = np.zeros(3) ; db[i] = h # set the field
            ds = np.zeros(3) ; ds[j] = h # set the field
            e00 = gete(-0.*db,-0.*ds) # --
            e01 = gete(-0.*db,ds) # -+
            e10 = gete(db,-0.*ds) # +-
            e11 = gete(db,ds) # +-
            dedb0 = (e10 - e00)/h
            dedb1 = (e11 - e01)/h
            o = (dedb1 - dedb0)/h # second derivative
            print(o,dedb0,e00,e01)
            out[i][j] = -2.*o # store
    return np.matrix(out) # return gtensor







