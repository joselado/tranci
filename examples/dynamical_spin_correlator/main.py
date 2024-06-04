import os ; import sys
cipath = os.getcwd() + "/../../" # address where the library is
sys.path.append(cipath+"/src")
import matplotlib.pyplot as plt

from tranci.atom import get_atom

Atom = get_atom(ne=6) # get the Atom
V = Atom.Operator["Coulomb"]
D = Atom.Operator["z2"] 

import numpy as np

H = 0.2*D + 4*V + 0.05*Atom.Operator["ls"]
M = Atom.get_manifolds(H)
A = Atom.SP_Operator["sy"]@Atom.SP_Operator["dz2"]
A = Atom.one2many(A) # to the many-body basis
print(M.get_gs_projected_eigenvalues(Atom.Operator["sz"]))
print(M.get_excitations())
(es,ds) = M.get_dynamical_correlator(A,es=np.linspace(-1,1.,100))

plt.plot(es,ds)

plt.show()

