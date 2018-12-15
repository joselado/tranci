
from __future__ import print_function
import numpy as np

def angular(x,y,z):
  xy = x*y - y*x
  xy = xy - 1j*z
  if np.abs(np.max(xy.todense()))>0.001: 
    print(x*y - y*x)
    raise



def zero(x,y):
  xy = x*y - y*x
  if np.abs(np.max(xy.todense()))>0.01: 
    print(x*y - y*x)
    raise












def check_all(at):
  """Check all the commutators"""
  # orbital
  angular(at.lx,at.ly,at.lz)
  angular(at.ly,at.lz,at.lx)
  angular(at.lz,at.lx,at.ly)
  # spin
  angular(at.sx,at.sy,at.sz)
  angular(at.sy,at.sz,at.sx)
  angular(at.sz,at.sx,at.sy)
  # ladder
  lp = at.lx + 1j*at.ly
  lm = at.lx - 1j*at.ly
#  print at.lz2
#  exit()
  angular(lp,lm,-2*1j*at.lz)
  # common
  zero(at.sx,at.lx)
  zero(at.sx,at.ly)
  zero(at.sx,at.lz)
  zero(at.sy,at.lx)
  zero(at.sy,at.ly)
  zero(at.sy,at.lz)
  zero(at.sz,at.lx)
  zero(at.sz,at.ly)
  zero(at.sz,at.lz)
  # coulomb 
  zero(at.lx,at.vc)
  zero(at.ly,at.vc)
  zero(at.lz,at.vc)
  zero(at.sx,at.vc)
  zero(at.sy,at.vc)
  zero(at.sz,at.vc)
