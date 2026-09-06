
from __future__ import print_function
# library with different hamiltonian
import scipy.sparse.linalg as lg
import scipy.linalg as dlg
import numpy as np

ntol = 6 # number of decimals to consider
ntol_ene = 6 # number of decimals to consider for the energies
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
  tri = p.trigonal
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
  if tri != 0.0: form += str(tri) + "\\left ( \\frac{l_x + l_y + l_z}{\\sqrt 3} \\right )^2  +"
  if z4 != 0.0: form += str(z4) + "l_z^4  +"
  if x2y2 != 0.0: form += str(x2y2) + "((l_xl_y)^2 + (l_yl_x)^2)  +"
  if soc != 0.0: form += str(soc) + "\\vec l \\cdot \\vec s  +"
#  if U != 0.0: form += str(U) + "V_{ijkl}c^\\dagger_i c^\\dagger_j c_k c_l  +"
  if U != 0.0: form += str(U) + "V_{e-e} +"  # U is a dimensionless multiplier
  if b[0] != 0.0: form += str(b[0]) + "(l_x+2s_x)  +"
  if b[1] != 0.0: form += str(b[1]) + "(l_y+2s_y)  +"
  if b[2] != 0.0: form += str(b[2]) + "(l_z+2s_z)  +"
  # exchange
  if j[0] != 0.0: form += str(j[0]) + "s_x  +"
  if j[1] != 0.0: form += str(j[1]) + "s_y  +"
  if j[2] != 0.0: form += str(j[2]) + "s_z  +"
  form = form.rstrip() # drop trailing whitespace
  if form.endswith("+"): form = form[:-1] # remove a dangling separator only
  form += "\\end{equation}\n" # end equation
  if n>-1: form += "Number of electrons in the d shell = "+str(n)+"\n\n"
  form += "Lower case $l,s$ denote single particle operators\n\n"
  form += "Upper case $L,S,J$ denote multi particle operators\n\n"
  return form




def build_hamiltonian(atom,p,gn=1./1836.):
  """Creates a simple Hamiltonian with D, E and soc parameters"""
  # get all the parameters
  # check if it is a dictionary
  # read the parameters the same way whether p is a dict or an object, so the
  # dictionary form does not silently miss O, z4 or b
  _missing = object()
  if isinstance(p,dict): get = lambda k: p.get(k,_missing)
  else: get = lambda k: getattr(p,k,_missing)
  def req(k): # a parameter the Hamiltonian cannot be built without
    v = get(k)
    if v is _missing: raise KeyError("build_hamiltonian: missing parameter '%s'"%k)
    return v
  def opt(k,default): # an optional parameter
    v = get(k)
    return default if v is _missing else v
  D = req("D")
  E = req("E")
  U = req("U")
  soc = req("soc")
  j = req("j")
  x2y2 = req("x2y2")
  trigonal = req("trigonal")
  O = opt("O",0.0)
  z4 = opt("z4",0.0)
  b = opt("b",[0.0,0.0,0.0])
  Uc = get("Uc")
  if Uc is _missing:
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
    theta = np.arccos(1./np.sqrt(3)) # theta angle
    phi = np.pi/4 # phi angle
    rot = lg.expm(1j*atom.lz*phi)@lg.expm(1j*atom.lx*theta) 
    h = h + trigonal*(rot@atom.z2@rot.conj().T) # trigonal field
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






class Lowest_States():
    """ Class for the lowest states"""
    has_degeneracies = False # if degeneracies have been calculated
    def __init__(self,h,atom=None):
        self.h = h # hamiltonian
        self.atom = atom # Atom object
        evals,evecs = eigenstates(h)
        self.e0 = float(np.min(evals)) # ground state energy before the shift
        evals = evals - np.min(evals)
        # Keep the exact eigenvalues. ntol/ntol_ene are degeneracy-grouping
        # tolerances, not a display precision: rounding here used to quantise
        # every reported energy at the grouping window (1e-3 eV from the GUI).
        self.evals = np.array(evals)
        self.evals_full = np.array(evals)
        self.evecs = evecs
    def get_representation(self,A,n=6):
        """Representation of a certain operator in a basis"""
        return get_representation(self.evecs[0:n],A)
    def get_gtensor(self,**kwargs):
        """Compute and return the g-tensor of a Kramers doublet ground state"""
        if self.atom is None:
            raise ValueError("get_gtensor needs the Atom object, build the "
                    "Lowest_States with atom=... or via Atom.get_manifolds")
        from .gtensor import get_gtensor
        self.gtensor = get_gtensor(self.atom,self.h,**kwargs)
        return self.gtensor
    def get_principal_g(self,**kwargs):
        """Return the principal g values and the magnetic axes"""
        if self.atom is None:
            raise ValueError("get_principal_g needs the Atom object, build the "
                    "Lowest_States with atom=... or via Atom.get_manifolds")
        from .gtensor import get_principal_g
        return get_principal_g(self.atom,self.h,**kwargs)
    def get_gs_degeneracy(self,tol=None,T=None):
      """Gets the degeneracy of the ground state manifold

      Returns (degeneracy, ground state energy). The degeneracy is the integer
      number of states within tol of the minimum, consistent with
      get_degeneracies/get_manifolds. It used to be a Boltzmann weight at a
      fixed absolute T=1e-4 eV, which is fractional for any manifold split on
      the sub-meV scale this code targets."""
      if tol is None: tol = globals()["tol"] if T is None else T # T: old name
      me = np.min(self.evals) # minimum energy
      de = np.abs(self.evals - me) # shift energy
      ngs = int(np.sum(de<tol)) # number of states within the tolerance
      return ngs,me
    def get_gs_multiplicity(self,**kwargs):
      """Get the ground state multiplicity"""
      return self.get_gs_degeneracy(**kwargs)[0]
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
      es = [es[i] - es[0] for i in range(len(es))] # return only en diff
      return es
    def get_gs_manifold(self):
      """Returns the vectors of the GS manifold"""
      self.gs_manifold = get_gs_manifold(self.evals,self.evecs)
      return self.gs_manifold
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
      return np.round(evals,6)
    def get_dynamical_correlator(self,A,B=None,**kwargs):
        if B is None: B = A
        from .dynamicstk import dynamics
        es,ds = 0,0
        gsm = self.get_gs_manifold() # ground state manifold
        e0 = self.e0 # unshifted ground state energy, computed once
        for wf0 in gsm:
            (ei,di) = dynamics.dynamical_correlator(self.h,
                    A,B,wf0=wf0,e0=e0,**kwargs)
            es = ei
            ds = ds + di
        return es,ds/len(gsm) # average over the degenerate ground states
    def get_correlation_entropy(self,wf):
        from . import entropy
        return entropy.correlation_entropy(self.atom,wf)


lowest_states = Lowest_States # alias for compatibility





def get_degeneracies(arr):
  """Get the degeneracies in an array"""
  me = np.min(arr) # minimum
  dg = 0
  arrrec = []
  for a in arr:
    if np.abs(a-me)<tol:
      dg += 1 # increase counter
    else:
      arrrec.append(a)
  pdg =(dg,me) # append degeneracy, with the exact manifold energy
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



def get_gs_manifold(evals,evecs,tol=None):
  """ Return a list with the GS manifold"""
  if tol is None: tol = globals()["tol"] # read at call time, as the siblings do
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
  # One sparse-times-dense product for all the vectors at once; the previous
  # implementation built three sparse matrices and densified per matrix element.
  W = np.array(wfs) # n x N matrix of wavefunctions
  return np.asmatrix(np.conjugate(W)@(A@W.T))



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


