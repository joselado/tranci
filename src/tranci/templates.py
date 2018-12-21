from . import hamiltonians
import numpy as np
from .numberformat import fform

def write_low_energy(ls,at,n=10,ghz=False):
  """Write energies and expectation values, input is lowest states"""
  waves = ls.evecs[0:n] # lowest states
  e0 = ls.evals[0:n] # lowest states
  # now write in a file
  h = 4.135*10**(-6) # in GHz
  fo = open("SPECTRUM.OUT","w")
  fo.write("# Energy [GHz], Energy [meV], Lx, Ly, Lz, Sx, Sy, Sz\n")
  for i in range(len(waves)): # loop over states
    fo.write(fform(e0[i]/h,n=6)+"    ") # write energy
    fo.write(fform(e0[i]*1000,n=6)+"    ") # write energy
    for op in [at.lx,at.ly,at.lz,at.sx,at.sy,at.sz]:
      a = hamiltonians.exp_val(waves[i],op)
      fo.write(fform(a)+"    ")
    fo.write("\n")
  fo.close()
  return np.genfromtxt("SPECTRUM.OUT") # return everything




def write_population(ls,at,n=10,ghz=False):
  """Write energies and expectation values, input is lowest states"""
  waves = ls.evecs[0:n] # lowest states
  e0 = ls.evals[0:n] # lowest states
  # now write in a file
  h = 4.135*10**(-6) # in GHz
  fo = open("POPULATION.OUT","w")
  fo.write("# Energy [GHz], Energy [meV], ")
  fo.write("up-2, up-1, up0, up+1, up+2, ")
  fo.write("dn-2, dn-1, dn0, dn+1, dn+2\n")
  for i in range(len(waves)): # loop over states
    fo.write(fform(e0[i]/h,n=6)+"    ") # write energy
    fo.write(fform(e0[i]*1000,n=6)+"    ") # write energy
    orbs = [at.um2,at.um1,at.u0,at.up1,at.up2]
    orbs += [at.dm2,at.dm1,at.d0,at.dp1,at.dp2]
    for op in orbs:
      a = hamiltonians.exp_val(waves[i],op)
      fo.write(fform(a)+"    ")
    fo.write("\n")
  fo.close()
  return np.genfromtxt("SPECTRUM.OUT") # return everything

