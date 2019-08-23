
from __future__ import print_function
# library with different hamiltonian
import scipy.sparse.linalg as lg
import scipy.linalg as dlg
import numpy as np

ntol = 12 # number of decimals to consider
ntol_ene = 12 # number of decimals to consider for the energies
tol = 10**(-ntol)
tol_ene = 10**(-ntol_ene)
scale_coulomb = 1.0 # constant to reproduce alejandro's results


def latex_DE(atom,p):
  """ Get the hamiltonian in latex form"""
  n = p.n
  D = p.D
  E = p.E
  U = p.U
  soc = p.soc
  O = p.O
  z4 = p.z4
  j = p.j
  b = p.b
  x2y2 = p.x2y2
#  form = "\\textit{Tranci version 0.3}\n\n\n"
  form = "\\nonstopmode\n\n\n" # ignore errors
  form += "\\section{Hamiltonian}\n" # name of the section
  form += "Hamiltonian of the atom\n"
  form += "\\begin{equation}\n" # begin equation
  form += "\\mathcal{H} = \\sum" # H
  if D != 0.0: form += str(D) + "l_z^2 + "
  if E != 0.0: form += str(E) + "(l_x^2 - l_y^2)  +"
  if O != 0.0: form += str(O) + "(l_x^4 + l_y^4 + l_z^4)  +"
  if z4 != 0.0: form += str(z4) + "l_z^4  +"
  if x2y2 != 0.0: form += str(x2y2) + "l_x^2l_y^2 + l_y^2l_x^2"
  if soc != 0.0: form += str(soc) + "\\vec l \\cdot \\vec s  +"
#  if U != 0.0: form += str(U) + "V_{ijkl}c^\\dagger_i c^\\dagger_j c_k c_l  +"
  if U != 0.0: form += str(U) + "V_{e-e} +"
  if b[0] != 0.0: form += str(b[0]) + "(l_x+2s_x)  +"
  if b[1] != 0.0: form += str(b[1]) + "(l_y+2s_y)  +"
  if b[2] != 0.0: form += str(b[2]) + "(l_z+2s_z)  +"
  # exchange
  if j[0] != 0.0: form += str(j[0]) + "s_x  +"
  if j[1] != 0.0: form += str(j[1]) + "s_y  +"
  if j[2] != 0.0: form += str(j[2]) + "s_z  +"
  form = form[:-1] # remove last character
  form += "\\end{equation}\n" # end equation
  if n>-1: form += "Number of electrons in the d shell = "+str(n)+"\n\n"
  form += "Lower case $l,s$ denote single particle operators\n\n"
  form += "Upper case $L,S,J$ denote multi particle operators\n\n"
  return form




def build_hamiltonian(atom,p,gn=1./1836.):
  """Creates a simple Hamiltonian with D, E and soc parameters"""
  # get all the parameters
  # check if it is a dictionary
  if isinstance(p,dict): # if it is a dictionary
    D = p["D"]
    E = p["E"]
    U = p["U"]
    soc = p["soc"]
    U = p["U"]
    j = p["j"]
    x2y2 = p["x2y2"]
    trigonal = p["trigonal"]
    try: cf = p["cf"]
    except: pass
    try: Uc = p["Uc"]
    except: 
      print("Using U as Uc")
      Uc = U
  else: # if it is a class
    D = p.D
    E = p.E
    U = p.U
    soc = p.soc
    O = p.O
    z4 = p.z4
    j = p.j
    b = p.b
    x2y2 = p.x2y2
    trigonal = p.trigonal
    try: cf = p.cf
    except: print("Not found total crystal field operator")
    try: Uc = p.Uc
    except: 
      print("Using U as Uc")
      Uc = U
  # build the hamiltonian
  h = soc*atom.ls + D*atom.z2 + E*(atom.x2-atom.y2)
  if Uc != 0.0:
    h = h + Uc*atom.vc*scale_coulomb 
    print("Added Coulomb using Coulomb integrals")
  if U != 0.0: # Add coulomb by symmetry
    try:
      h = h - U*atom.dl2/3. # maximize L2
      h = h - U*atom.ds2 # maximize S2
      print("Added Coulomb by symmetry")
    except: print("No Coulomb added by symmetry")
#  h = soc*atom.ls + D*atom.z2 + E*(atom.x2-atom.y2)
  h = h + O*(atom.x4 + atom.y4 + atom.z4) # octahedral field
  h = h + z4*atom.z4 # octahedral field
  h = h + x2y2*atom.x2y2 # square field
  if trigonal != 0.0:
    rot = lg.expm(1j*atom.lz*np.pi/4)*lg.expm(1j*atom.lx*np.pi/4) 
    h = h + trigonal*(rot*atom.z2*rot.H) # trigonal field
  h = h + j[0]*(atom.sx)  # Zeeman x
  h = h + j[1]*(atom.sy)  # Zeeman y
  h = h + j[2]*(atom.sz)  # Zeeman z
  h = h + b[0]*(atom.lx + 2*atom.sx)  # Zeeman x
  h = h + b[1]*(atom.ly + 2*atom.sy)  # Zeeman y
  h = h + b[2]*(atom.lz + 2*atom.sz)  # Zeeman z
  try:
    h = h + atom.cf
    print("Added total crystal field")
  except: print("Total CF not added")
  # hyperfine coupling
  if atom.has_nucleus: 
    h = h + p.hyper*atom.si # hyperfine coupling 
    h = h + gn*(p.b[0]*atom.ix + p.b[1]*atom.iy + p.b[2]*atom.iz) # Zeeman coupling
  return h




def eigenstates(h,n=20,maxiter=None):
#  (evals,evecs) = lg.eigsh(h,k=n,which="SA",maxiter=10000)
  h2 = h.todense()*1000.0 # to meV
  (evals,evecs) = dlg.eigh(h2)
  evals = evals/1000.0 # to eV
  evecs = evecs.transpose() # transpose eigenvectors
#  evecs = [v for (e,v) in sorted(zip(evals,evecs))] # sort eigenvectors
#  evals = sorted(evals) # sort eigenvalues
  return (evals,evecs)



def eigenvalues(h,n=20):
  """ Return the eigenvalues"""
  return eigenstates(h,n=n)[0]






class lowest_states():
  """ Class for the lowest states"""
  has_degeneracies = False # if degeneracies have been calculated
  def __init__(self,h,n=10,atom=None):
    self.h = h # hamiltonian
    self.atom = atom # Atom object
    evals,evecs = eigenstates(h,n=10)
    evals = evals - np.min(evals)
    self.evals = np.array([round(e,ntol) for e in evals]) # round values
    self.evals_full = np.array([round(e,ntol_ene) for e in evals]) # round values
    self.evecs = evecs
  def get_gtensor(self):
    """Compute the gtensor"""
    if self.atom is None: raise
    from .gtensor import get_gtensor
    self.gtensor = get_gtensor(self.atom,self.h)
  def get_gs_degeneracy(self):
    """Gets the degeneracy of each manifold"""
    me = np.min(self.evals) # minimum energy
    de = np.abs(self.evals - me) # shift energy
    ngs = len(de[de<=tol]) # number of states within an interval
    return ngs,me
  def get_gs_multiplicity(self):
    """Get the ground state multiplicity"""
    return self.get_gs_degeneracy()[0]
  def project_operator(self,m):
    """ Gets the proyection of an operator of the low energy states"""
  def get_degeneracies(self):
    """ Gets the degeneracies of the states diagonalized"""
    if not self.has_degeneracies: # if not calculated yet
      self.degeneracies = get_degeneracies(self.evals)
      self.has_degeneracies = True
    return self.degeneracies
  def get_multiplicities(self):
    """ Get multiplicity of the manifolds"""
    dgs = self.get_degeneracies()
    return [d[0] for d in dgs] # return only the degeneracies
  def get_excitations(self):
    """Gets the energies of the excited states"""  
    dgs = self.get_degeneracies()
    es = [d[1] for d in dgs] # return only the eigenvalues
    es = [round(es[i] - es[0],ntol) for i in range(len(es))] # return only en diff
    return es
  def get_gs_manifold(self):
    """Returns the vectors of the GS manifold"""
    self.gs_manifold = get_gs_manifold(self.evals,self.evecs)
  def get_manifolds(self):
    """ Returns a list with the different manifolds"""
    self.manifolds = get_manifolds(self.evals,self.evecs) # store in object
    return self.manifolds # return the manifolds
  def disentangle_manifolds(self,a):
    """ Disentangle the manifolds according to an operator"""
    self.get_manifolds() # get the manifolds
    mani = [disentangle_manifold(wfl,a) for wfl in self.manifolds]
    self.manifolds = mani # put new manifolds
    self.gs_manifold = mani[0] # put new manifold
    waves = [] # empty list
    for m in mani:
        for w in m: waves.append(w)
    self.evecs = np.array(waves) # store disentangled waves
  def disentangle_gs_manifold(self,a):
    """ Disentangle the manifolds according to an operator"""
    self.get_gs_manifold() # get the manifolds
    self.gs_manifold = disentangle_manifold(self.gs_manifold,a) 
  def get_gs_projected_eigenvalues(self,A):
    """ Get the projected eigenvalues of a certain operator"""
    self.get_gs_manifold() # get the manifold
    evals = get_projected_eigenvalues(self.gs_manifold,A)  # diagonalize
    return evals








def get_degeneracies(arr):
  """Get the degeneracies in an array"""
  me = min(arr) # minimum
  dg = 0
  arrrec = []
  for a in arr:
    if abs(a-me)<tol:
      dg += 1 # increase counter
    else:
      arrrec.append(a)
  pdg =(dg,round(me,ntol)) # append degeneracy
  if len(arrrec)>0: 
    return [pdg] + get_degeneracies(arrrec) # if still numbers, iterate
  else: 
    return [pdg]  # if the remaining list is empty, return 


def get_manifolds(evals,evecs):
  """ Return a list with the different manifolds, splitted
  by energy """
  me = min(evals) # minimum
  wfm = [] # list for the wavefunctions in this manifold
  evalsrec = [] # list for eigenvalues left
  evecsrec = [] # list for eigenfunctions left
  for (a,v) in zip(evals,evecs):
    if abs(a-me)<tol:
      wfm.append(v) # append wavefunction
    else:
      evalsrec.append(a) # store eigval left
      evecsrec.append(v) # eigfun left
  if len(evalsrec)>0: 
    return [wfm] + get_manifolds(evalsrec,evecsrec) # if still wfs, iterate
  else: 
    return [wfm] # return this manifold



def get_gs_manifold(evals,evecs,tol=tol):
  """ Return a list with the GS manifold"""
  me = min(evals) # minimum
  wfm = [] # list for the wavefunctions in this manifold
  evalsrec = [] # list for eigenvalues left
  evecsrec = [] # list for eigenfunctions left
  for (a,v) in zip(evals,evecs):
    if abs(a-me)<tol:
      wfm.append(v) # append wavefunction
  return wfm # return this manifold








def get_representation(wfs,A):
  """Gets the matrix representation of a certain operator"""
  n = len(wfs) # number of eigenfunctions
  ma = np.matrix([[0.0j for i in range(n)] for j in range(n)]) # representation of A
  from scipy.sparse import csc_matrix as csc
  sa = csc(A) # sparse matrix
  for i in range(n):
    vi = csc(np.conjugate(wfs[i])) # first wavefunction
    for j in range(n):
      vj = csc(wfs[j]).transpose() # first wavefunction
      data = (vi*sa*vj).todense()[0,0]
      ma[i,j] = data
  return ma



def exp_val(wf,A):
  """Gets the matrix representation of a certain operator"""
  from scipy.sparse import csc_matrix as csc
  sa = csc(A) # sparse matrix
  vi = csc(np.conjugate(wf)) # first wavefunction
  vj = csc(wf).transpose() # first wavefunction
  data = (vi*sa*vj).todense()[0,0].real
  return data





def get_projected_eigenvalues(wfs,A):
  """Get eigenvalues of a proyected operator"""
  ma = get_representation(wfs,A)
  evals = dlg.eigvalsh(ma)
  return evals  # return eigenvalues





def disentangle_manifold(wfs,A):
  """ Disentangles the wavefunctions of a degenerate manifold
  by expressing them in terms of eigenvalues of an input operator"""
  ma = get_representation(wfs,A) # get the matrix form of the operator
  wfsout = [] # empty list
  evals,evecs = dlg.eigh(ma) # diagonalize
  evecs = evecs.transpose() # transpose eigenvectors
  for v in evecs: # loop over eigenvectors
    wf = wfs[0]*0.0
    for (i,iv) in zip(range(len(v)),v): # loop over components
      wf += iv*wfs[i] # add contribution
    wfsout.append(wf.copy()) # store wavefunction
  return wfsout


