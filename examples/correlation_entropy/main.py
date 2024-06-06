import os ; import sys
cipath = os.getcwd() + "/../../" # address where the library is
sys.path.append(cipath+"/src")
import matplotlib.pyplot as plt

from tranci.atom import get_atom

Atom = get_atom(ne=3) # get the Atom
V = Atom.Operator["Coulomb"]
D = Atom.Operator["z2"] 
Sz = Atom.Operator["sz"] 
LS = Atom.Operator["ls"] 

H = 4*V + D + 0.01*Sz + 0.05*LS
M = Atom.get_manifolds(H)
S = M.get_correlation_entropy(M.get_gs_manifold()[0])
print(S)
