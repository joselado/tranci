# enlarge to consider hyperfine coupling
from copy import deepcopy
from scipy.sparse import csc_matrix,bmat
import numpy as np

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
    iden = csc_matrix(np.identity(norig,dtype=np.complex)) # identity operator
    iz = bmat([[iden,None],[None,-iden]])/2.
    ix = bmat([[None,iden],[iden,None]])/2.
    iy = bmat([[None,-1j*iden],[1j*iden,None]])/2.
  else: raise # raise error
  for key in at.terms: # loop over 
    at.terms[key] = enlarge(at.terms[key],n=ns) # enlarge Hamiltonian
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
    



  

    






