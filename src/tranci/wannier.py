import angular
from scipy.sparse import csc_matrix as csc
from scipy.sparse import coo_matrix
from scipy.sparse import bmat
from scipy.optimize import minimize
import numpy as np


def get_cf(hname="cf.in",oname="worbitals.in"):
  """Read the matrix in the correct order"""
  # first read the matrix
  mm = np.genfromtxt(hname).transpose() # matrix of the crystal field
  norb = int(max(mm[0])) # number of orbitals
  # total dimension
  m = np.matrix([[0.0j for i in range(norb)] for j in range(norb)]) 
  for i in range(len(mm[0])):
    m[int(mm[0][i]-1),int(mm[1][i]-1)] = mm[2][i] + 1j*mm[3][i]
    
  # now sort the matrix   
  inorbs = open(oname).readlines() # read the lines
  inorbs = [d.split()[0] for d in inorbs] # clean
  dorbs = ["dz2","dxz","dyz","dxy","dx2-y2"] # this the right order

  # Now read the rest of the orbitals
  nf = norb - 5 # number of fluctuating
  print(nf,"fluctuating orbitals") 
  forbs = [] # empty list
  for i in inorbs: # loop over all orbitals
    if i not in dorbs: forbs.append(i) # store fluctuating orbital
  torbs = dorbs + forbs # total orbitals, first d
  print("List of orbitals",torbs)
  inds = dict() # create dictionary
  for d in torbs: # loop over orbitals
    for i in range(len(inorbs)): # loop over input orbitals
      if d in inorbs[i]: inds[d] = i # store this index
  # now read the matrix in the right order
  mout = np.matrix([[0.0j for i in range(len(torbs))] for j in range(len(torbs))])


  for i in range(len(torbs)):
    for j in range(len(torbs)):
      ii,jj = inds[torbs[i]], inds[torbs[j]] # get the right indexes
      mout[i,j] = m[ii,jj] # get the right number

  Rd = csc(ylm2xyz_l2()) # get the rotation matrix
  Ri = csc(np.identity(nf)) # identity matrix
  R = bmat([[Rd,None],[None,Ri]]).todense()
  cf = R*mout*R.H # return the matrix in spherical harmonics
  # now write the matrix in a file so that the fluctuating program can read it
  sm = coo_matrix(cf) # in sparse form
  fw = open("hopping.in","w")
  fw.write(str(len(sm.data))+"\n")
  data = sm.data
  for (i,j,d) in zip(sm.row,sm.col,data):
    fw.write(str(i)+"   ")
    fw.write(str(j)+"   ")
    fw.write(str(d.real)+"   ")
    fw.write(str(d.imag)+"\n")
  fw.close()
  print("Written crystal field in hopping.in")



def wannier2cf(m,terms=None,cftype="all"):
  """Transform the Wannier basis into a crystal field Ylm one"""
  if terms==None:
    (terms,names) = get_operators(cftype) # get the operators
  def fun(x):
    dif = m*0. # initialice
    for i in range(len(x)):
      dif += x[i]*terms[i] # add term
      r = m-dif # difference
      r = r*r.H # all eigenvalues are now positive
    return r.trace()[0,0] # return the deviation 
  x0 = np.random.random(len(terms))-0.5 # initialice
  result = minimize(fun,x0) # minimize the function
  outcf = dict() # create a dictionary
  for (t,n,ix) in zip(terms,names,result.x):
    outcf[n] = (ix,t)
    print(n,ix)
  print("Error",fun(result.x))
  return outcf # return the dictionary
  # fit to  



def get_operators(cftype):
  """Return the operators of the CF"""
  (lx,ly,lz) = angular.angular(l=2) # angular terms
  terms = [] # empty list
  if cftype=="all":
    terms = [lx**4,ly**4,lz**4] # quartic
    terms += [ly**2*lz**2,lx**2*lz**2,lx**2*ly**2] # quartic C4
    names = ["x4+y4+z4","y2z2","x2z2","x2y2"]
  elif cftype=="C4":
    terms = [np.matrix(np.identity(5))] # square
    terms += [lz**2] # square
    terms += [lz**4] # quartic
    terms += [lx**4 + ly**4 + lz**4] # octahedral
    names = ["I","z2","z4","x4+y4+z4"] # names
  else: raise
  return terms,names



def ylm2xyz_l2():
  """Return the matrix that converts the cartesian into spherical harmonics"""
  m = np.matrix([[0.0j for i in range(5)] for j in range(5)])
  s2 = np.sqrt(2.)
  m[2,0] = 1. # dz2
  m[1,1] = 1./s2 # dxz
  m[3,1] = -1./s2 # dxz
  m[1,2] = 1j/s2 # dyz
  m[3,2] = 1j/s2 # dyz
  m[0,3] = 1j/s2 # dxy
  m[4,3] = -1j/s2 # dxy
  m[0,4] = 1./s2 # dx2y2
  m[4,4] = 1./s2 # dx2y2
  return m # return change of bassi matrix

