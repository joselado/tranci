import numpy as np
import cmath
from .read import read_matrix
from .read import read_basis
from .numberformat import zform

guess_number = True
tol = 0.01


class CIatom():
  """Class to work with a CI atom"""
  has_nucleus = False # does not have nucleus
  def read(self,path=""):
    """ Read all the matrices of the hamiltonian"""
    try:  self.cf = read_matrix(path+"hopping.op")   # read Sx
    except: print("hopping.op not found")
    try:  self.ds2 = read_matrix(path+"ds2.op")   # read Sx
    except: print("ds2.op not found")
    try:  self.dl2 = read_matrix(path+"dl2.op")   # read Sx
    except: print("dl2.op not found")
    self.sx = read_matrix(path+"sx.op")   # read Sx
    self.sy = read_matrix(path+"sy.op")   # read Sy
    self.sz = read_matrix(path+"sz.op")   # read Sz
    self.s2 = self.sx*self.sx + self.sy*self.sy + self.sz*self.sz
    self.lx = read_matrix(path+"lx.op")   # read Lx
    self.ly = read_matrix(path+"ly.op")   # read Ly
    self.lz = read_matrix(path+"lz.op")   # read Lz
    self.l2 = self.lx*self.lx + self.ly*self.ly + self.lz*self.lz
    self.z2 = read_matrix(path+"z2.op")   # read Lz
    self.x2 = read_matrix(path+"x2.op")   # read Lz
    self.y2 = read_matrix(path+"y2.op")   # read Lz
    self.x2y2 = read_matrix(path+"x2y2.op")   # read Lz
    self.z4 = read_matrix(path+"z4.op")   # read Lz
    self.x4 = read_matrix(path+"x4.op")   # read Lz
    self.y4 = read_matrix(path+"y4.op")   # read Lz
    self.jx = read_matrix(path+"jx.op")   # read Jx
    self.jy = read_matrix(path+"jy.op")   # read Jy
    self.jz = read_matrix(path+"jz.op")   # read Jz
    self.j2 = self.jx*self.jx + self.jy*self.jy + self.jz*self.jz
    self.vc = read_matrix(path+"vc.op")   # read coulomb term
    self.ls = read_matrix(path+"ls.op")   # read coulomb term
    # read all the projection operators
    self.um2 = read_matrix(path+"0.op")
    self.um1 = read_matrix(path+"1.op")
    self.u0 = read_matrix(path+"2.op")
    self.up1 = read_matrix(path+"3.op")
    self.up2 = read_matrix(path+"4.op")
    self.dm2 = read_matrix(path+"5.op")
    self.dm1 = read_matrix(path+"6.op")
    self.d0 = read_matrix(path+"7.op")
    self.dp1 = read_matrix(path+"8.op")
    self.dp2 = read_matrix(path+"9.op")
    # add a dictionary
    terms = dict()
    terms["sx"] = self.sx
    terms["sy"] = self.sy
    terms["sz"] = self.sz
    terms["lx"] = self.lx
    terms["ly"] = self.ly
    terms["lz"] = self.lz
    terms["jx"] = self.jx
    terms["jy"] = self.jy
    terms["jz"] = self.jz
    terms["s2"] = self.s2
    terms["l2"] = self.l2
    terms["j2"] = self.j2
    terms["x2"] = self.x2
    terms["y2"] = self.y2
    terms["z2"] = self.z2
    terms["x4"] = self.x4
    terms["y4"] = self.y4
    terms["z4"] = self.z4
    terms["x2y2"] = self.x2y2
    terms["vc"] = self.vc
    terms["ls"] = self.ls
    try: terms["cf"] = self.cf
    except: pass
    self.terms = terms
  def update(self):
      self.sx  = self.terms["sx"]
      self.sy  = self.terms["sy"]
      self.sz  = self.terms["sz"]
      self.lx  = self.terms["lx"]
      self.ly  = self.terms["ly"]
      self.lz  = self.terms["lz"]
      self.jx  = self.terms["jx"]
      self.jy  = self.terms["jy"]
      self.jz  = self.terms["jz"]
      self.s2  = self.terms["s2"]
      self.l2  = self.terms["l2"]
      self.j2  = self.terms["j2"]
      self.x2  = self.terms["x2"]
      self.y2  = self.terms["y2"]
      self.z2  = self.terms["z2"]
      self.x4  = self.terms["x4"]
      self.y4  = self.terms["y4"]
      self.z4  = self.terms["z4"]
      self.x2y2 = self.terms["x2y2"]
      self.vc  = self.terms["vc"]
      self.ls  = self.terms["ls"]
  def get_basis(self,path=""):
    self.basis = read_basis("basis.out",path) # read the basis
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




def format_component(x):
  """Returns a string with the number correctly formatted"""
  if np.abs((np.abs(x)-1.))<tol: return ""
  else: return zform(x)




def get_atom(ne=None):
  """Get an atom, by the number of electrons"""
  import os  
  if ne is None:
    path = os.getcwd()+"/" #
  else:
    path = os.environ["TRANCIROOT"]+"/cilib/"+str(ne)+"/" #
  print("Reading from",path)
  at = CIatom() # create the CI object
  at.read(path=path) # read all the matrices
  at.get_basis(path=path) # read the basis from file
  return at # return atom


