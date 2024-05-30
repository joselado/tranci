import numpy as np
import cmath

guess_number = True
tol = 0.01


class CIatom():
  """Class to work with a CI atom"""
  def read(self):
    """ Read all the matrices of the hamiltonian"""
    from read import read_matrix
    self.sx = read_matrix("sx.op")   # read Sx
    self.sy = read_matrix("sy.op")   # read Sy
    self.sz = read_matrix("sz.op")   # read Sz
    self.s2 = self.sx*self.sx + self.sy*self.sy + self.sz*self.sz
    self.lx = read_matrix("lx.op")   # read Lx
    self.ly = read_matrix("ly.op")   # read Ly
    self.lz = read_matrix("lz.op")   # read Lz
    self.l2 = self.lx*self.lx + self.ly*self.ly + self.lz*self.lz
    self.z2 = read_matrix("z2.op")   # read Lz
    self.x2 = read_matrix("x2.op")   # read Lz
    self.y2 = read_matrix("y2.op")   # read Lz
    self.z4 = read_matrix("z4.op")   # read Lz
    self.x4 = read_matrix("x4.op")   # read Lz
    self.y4 = read_matrix("y4.op")   # read Lz
    self.jx = read_matrix("jx.op")   # read Jx
    self.jy = read_matrix("jy.op")   # read Jy
    self.jz = read_matrix("jz.op")   # read Jz
    self.j2 = self.jx*self.jx + self.jy*self.jy + self.jz*self.jz
    self.vc = read_matrix("vc.op")   # read coulomb term
    self.ls = read_matrix("ls.op")   # read coulomb term
  def get_basis(self):
    from read import read_basis
    self.basis = read_basis("basis.out") # read the basis
  def get_latex_wavefunction(self,wf):
    """ Outputs a wavefunction in latex format"""
    tol = 0.01
    strwf = [] # string for the wavefunction
    norm = []
    for (iwf,b) in zip(wf,self.basis):
      if np.abs(iwf*iwf)>tol:
        strwf.append(format_component(iwf)+"  &  "+b.get_latex())
        norm.append(-np.abs(iwf))
    strwf = [y for (x,y) in sorted(zip(norm,strwf))]
    return strwf # return a list of contributions



from format import zform

def format_component(x):
  """Returns a string with the number correctly formatted"""
  if np.abs((np.abs(x)-1.))<tol: return ""
  else: return zform(x)

