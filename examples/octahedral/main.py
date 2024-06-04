import os ; import sys
cipath = os.getcwd() + "/../../" # address where the library is
sys.path.append(cipath+"/src")
import matplotlib.pyplot as plt

from tranci.atom import get_atom

Atom = get_atom(ne=1) # get the Atom
lx = Atom.Operator["lx"]
ly = Atom.Operator["ly"]
lz = Atom.Operator["lz"]
V = Atom.Operator["Coulomb"]
H = lx@lx@lx@lx + ly@ly@ly@ly +lz@lz@lz@lz + 0.1*ly@ly + 4*V

M = Atom.get_manifolds(H)
print("Degeneracy",M.get_gs_degeneracy()[0])
#print("Spin",M.get_gs_projected_eigenvalues(Atom.Operator["sz"]))
for o in ["dz2","dxz","dyz","dx2y2","dxy"]:
#    print(o)
#    print(Atom.Operator[o])
    print(o,M.get_gs_projected_eigenvalues(Atom.Operator[o]))

