# routines to compute the gtensor

def get_gtensor(atom,h0):
    """Compute the gtensor according to the definition g = d2E/dBdS"""
    def geth(h0,db,ds): # get Hamiltonian for infinitesimal db and ds
        out = h0 # first term
        out = out + db[0]*(atom.lx + 2.*ds[0]*atom.sx) # first component
        out = out + db[1]*(atom.ly + 2.*ds[1]*atom.sy) # second component
        out = out + db[2]*(atom.lz + 2.*ds[2]*atom.sz) # third component
        return out # return Hamiltonian





