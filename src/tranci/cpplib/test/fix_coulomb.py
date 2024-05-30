tol = 0.01

import numpy as np

def fix_coulomb(atom):
  """ Guarantees that the coulomb term has rotational symmetry"""


def com(a,b):
  """ Returns the commutatow"""
  return a*b - b*a



def solve_symmetrizazion(oin,osy):
  """ Symmetrizes a certain operator with respect another one"""
  oout = oin + com(oin,osy) # selfconsistent equation
  error = np.max(np.abs((oout-oin).todense()))
  if error<tol: return oin
  else: 
    print "Error=",error
    solve_symmetrizazion((oout+oin)/2.0,osy)
  

