import numpy as np


def angular(l=2):
  """Calculate single particle angular momentum"""
  nm = 2*l + 1 # number of components
  zero = np.matrix([[0.0j for i in range(nm)] for j in range(nm)])
  # initialize matrices
  lz = zero.copy()
  lx = zero.copy()
  ly = zero.copy()
  lm = zero.copy()
  lp = zero.copy()
  # create l+ and l- and lz
  for m in range(-l,l): # loop over m components
    val = np.sqrt((l-m)*(l+m+1)) # value of the cupling
    im = m + l
    lp[(im+1),im] = val # up channel
  for m in range(-l,l+1): # loop over m components
    im = m + l
    lz[im,im] = m # value of lz, up channel
  lm = lp.H # adjoint
  lx = (lp + lm) /2.
  ly = -1j*(lp - lm) /2.
  return (lx,ly,lz)
