import os ; import sys
cipath = os.getcwd() + "/../../" # address where the library is
sys.path.append(cipath+"/src")
import matplotlib.pyplot as plt

from tranci.atom import get_atom

Atom = get_atom(ne=2) # get the Atom
lz = Atom.Operator["lz"]
V = Atom.Operator["Coulomb"]
H = lz@lz + 4*V

M = Atom.get_manifolds(H)
print("Degeneracy",M.get_gs_degeneracy()[0])
print("Spin",M.get_gs_projected_eigenvalues(Atom.Operator["sz"]))

