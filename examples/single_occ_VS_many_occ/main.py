import os ; import sys
cipath = os.getcwd() + "/../../" # address where the library is
sys.path.append(cipath+"/src")
import matplotlib.pyplot as plt

from tranci.atom import get_atom

# see how the occupation of a many-body atom
# changes as a function of the occupation
# of a single particle one

Atom = get_atom(ne=2) # get the Atom
Atom_single = get_atom(ne=1) # get the Atom
H0_MB = Atom.Operator["Coulomb"]
x4 = Atom_single.Operator["x4"] 
y4 = Atom_single.Operator["y4"] 
z4 = Atom_single.Operator["z4"] 
H0_single = 0.1*(x4 + y4 + z4) + 0.05*Atom_single.Operator["ls"]
V_single = Atom_single.Operator["z2"]  # perturbation

import numpy as np

CFs = np.linspace(0.1,2.,20)
io = 0 # counter

fig = plt.figure(figsize=[14,5])

plt.subplot(1,2,1)


# many body Hamiltonian occupations

plt.title("Many body")
io = 0
for o in ["dz2","dxz","dyz","dx2y2","dxy"]:
    io += 1
    ns = []
    for C in CFs:
        H = Atom.one2many(H0_single + C*V_single) + H0_MB 
        M = Atom.get_manifolds(H)
        ni = M.get_gs_projected_eigenvalues(Atom.Operator[o]) # projections
        ns.append(np.mean(ni)) # average occupation
    plt.scatter(CFs,ns,label=o,s=20*(6-io))
    plt.legend()
    plt.xlabel("Crystal field")
    plt.ylabel("Occupation")

# single particle occupations
plt.subplot(1,2,2)

plt.title("Single particle")
io = 0
for o in ["dz2","dxz","dyz","dx2y2","dxy"]:
    io += 1
    ns = []
    for C in CFs:
        H = Atom_single.one2many(H0_single + C*V_single)
        M = Atom_single.get_manifolds(H)
        ni = M.get_gs_projected_eigenvalues(Atom_single.Operator[o]) # proj
        ns.append(np.mean(ni)) # average occupation
    plt.scatter(CFs,ns,label=o,s=20*(6-io))
    plt.legend()
    plt.xlabel("Crystal field")
    plt.ylabel("Occupation")


plt.tight_layout()
plt.show()
