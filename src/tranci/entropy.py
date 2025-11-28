import numpy as np
from .dynamicstk import algebra

def correlation_entropy(Atom,wf):
    """Compute the correlation entropy of a wavefunction"""
    norb = 10 # number of orbitals
    dm = np.zeros((norb,norb),dtype=np.complex128)
    for i in range(norb):
        for j in range(i,norb): # loop
            m0 = np.zeros((norb,norb),dtype=np.complex128)
            m0[i,j] = 1.0 # this element
            m = Atom.one2many(m0) # many body
            d0 = algebra.braket_wAw(wf,m) # matrix element
            dm[i,j] = d0
            dm[j,i] = np.conjugate(d0)
    es = algebra.eigvalsh(dm) # eigenvalues
    es = es[es>1e-7] # only positive
    return -np.sum(es*np.log(es))



