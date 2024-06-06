import os ; import sys
cipath = os.getcwd() + "/../../" # address where the library is
sys.path.append(cipath+"/src")
import matplotlib.pyplot as plt

from tranci.atom import get_atom
import numpy as np

ne = 2
# generate an atom class for the system you want to compute
Atom = get_atom(ne=ne)
# total Hamiltonian
V = Atom.Operator["Coulomb"]
CF = Atom.Operator["z2"] #+ 0.2*(Atom.Operator["x2"] - Atom.Operator["y2"])
LS = Atom.Operator["ls"]
H0 = 4*V + 0.3*CF  # original Hamiltonian
H = H0 + 0.1*LS # original Hamiltonian
# degeneracy without magnetic field
M0 = Atom.get_manifolds(H0)
M = Atom.get_manifolds(H)

deg = M0.get_gs_multiplicity()
print("Unperturbed multiplicity",deg)
deg1 = M.get_gs_multiplicity()
print("New multiplicity",deg1)
deg = int(np.round(deg,1)) # degeneracy (as integer)

print("Sz",M0.get_gs_projected_eigenvalues(Atom.Operator["sz"]))
print("Lz",M0.get_gs_projected_eigenvalues(Atom.Operator["lz"]))

from tranci import effectivehamiltonian
#deg = 4
print("Hamiltonian projected on the ",deg," lowest states")
text = effectivehamiltonian.effective_spin_hamiltonian(M,H=H,n=deg,
        operators=["sx","sy","sz"])
#        operators=["sx","sy","sz","lx","ly","lz"])

print(text)

#from IPython.display import display, Math, Latex
#Latex(text)
