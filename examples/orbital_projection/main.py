import os ; import sys
cipath = os.getcwd() + "/../../" # address where the library is
sys.path.append(cipath+"/src")
import matplotlib.pyplot as plt

from tranci.atom import get_atom

Atom = get_atom(ne=3) # get the Atom
V = Atom.Operator["Coulomb"]
D = Atom.Operator["z2"] 

import numpy as np

CFs = np.linspace(0.1,2.,30)
io = 0 # counter

fig = plt.figure(figsize=[14,5])

plt.subplot(1,3,1)

for o in ["dz2","dxz","dyz","dx2y2","dxy"]:
    io += 1
    ns = []
    for C in CFs:
        H = C*D + 4*V
        M = Atom.get_manifolds(H)
        ni = M.get_gs_projected_eigenvalues(Atom.Operator[o]) # projections
        ns.append(np.mean(ni)) # average occupation
#    plt.subplot(3,2,io)
    plt.scatter(CFs,ns,label=o)
    plt.legend()
    plt.xlabel("Crystal field")
    plt.ylabel("Occupation")

# degeneracies

plt.subplot(1,3,3)
for C in CFs:
    H = C*D + 4*V
    M = Atom.get_manifolds(H)
    szi = M.get_gs_projected_eigenvalues(Atom.Operator["sz"])
    plt.scatter(szi*0. + C,szi,c="red")

plt.xlabel("Crystal field")
plt.ylabel("Sz")


# degeneracies
plt.subplot(1,3,2)
degs = []
for C in CFs:
    H = C*D + 4*V
    deg = Atom.get_manifolds(H).get_gs_degeneracy()[0]
    degs.append(deg)
plt.scatter(CFs,degs)

plt.xlabel("Crystal field")
plt.ylabel("GS degeneracy")



plt.tight_layout()
plt.show()
