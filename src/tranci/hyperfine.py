# enlarge to consider hyperfine coupling
from copy import deepcopy
from scipy.sparse import csc_matrix,bmat
import numpy as np

def enlarge_basis(basis,ns):
  """Return the product basis of the electronic basis and the nuclear spin

  enlarge() builds block(m,m), so state i of block b is |basis[i]> x |nuclear b>.
  Without this the operators grow while atom.basis does not, and
  get_latex_wavefunction silently drops half of every wavefunction."""
  if ns!=2: raise ValueError("Only spin 1/2 nuclei are implemented")
  labels = ["I\\uparrow","I\\downarrow"] # nuclear states, in block order
  out = [] # new basis
  for ib in range(ns): # loop over nuclear states
    for b in basis: # loop over electronic states
      bo = deepcopy(b) # copy the electronic state
      nuc = [0]*ns # occupation of the nuclear levels
      nuc[ib] = 1 # this nuclear state is occupied
      bo.v = np.concatenate([np.array(b.v),np.array(nuc)]) # extended vector
      bo.orb = list(b.orb) + labels # extended orbital labels
      bo.num2occ(bo.v) # recompute the occupations
      out.append(bo) # store
  return out


def add_nucleus(atin,s=0.5):
  at = deepcopy(atin) # copy object
  norig = atin.sx.shape[0] # original shape
  def enlarge(m,n=2):
    mo = [[None for j in range(n)] for i in range(n)]
    for i in range(n): mo[i][i] = m
    mo = bmat(mo) # return matrix
    return mo
  if s==0.5: 
    ns = 2 # additional size, nuclear spin
    iden = csc_matrix(np.identity(norig,dtype=np.complex128)) # identity operator
    iz = bmat([[iden,None],[None,-iden]])/2.
    ix = bmat([[None,iden],[iden,None]])/2.
    iy = bmat([[None,-1j*iden],[1j*iden,None]])/2.
  else: raise # raise error
  # rebind the entries of a fresh dictionary rather than mutating the one the
  # input atom (and its SP_Operator) may still be sharing
  at.terms = dict(at.terms)
  for key in at.terms: # loop over 
    at.terms[key] = enlarge(at.terms[key],n=ns) # enlarge Hamiltonian
  at.Operator = at.terms # keep the documented alias pointing at the new dict
  at.basis = enlarge_basis(atin.basis,ns) # the basis must grow with the operators
  at.update() # update attributes
  # add the new atributes
  at.ix = ix 
  at.iy = iy 
  at.iz = iz
  at.terms["ix"] = ix
  at.terms["iy"] = iy
  at.terms["iz"] = iz
  # add a new attribute
  si = at.terms["sx"]*ix + at.terms["sy"]*iy + at.terms["sz"]*iz
  at.si = si  # add hyper
  at.terms["si"] = si # add hyperfine coupling
  at.has_nucleus = True # has nucleus
  return at
    



  

    






