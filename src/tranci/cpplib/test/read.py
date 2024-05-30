import numpy as np

def read_matrix(namefile):
  """ Function which reads a sparse matrix from file"""
  d = open(namefile,"r").readlines()[0] # first line
  d = int(d.split("=")[1]) # dimension of the matrix
  from scipy.sparse import csc_matrix
  try:
    m = np.genfromtxt(namefile) # read the file
    m = m.transpose() # transpose matrix
    row = [int(i) for i in m[0]] # row
    col = [int(i) for i in m[1]] # column
    data = [i+1j*j for (i,j) in zip(m[2],m[3])] # data
    mout = csc_matrix((data,(row,col)),shape=(d,d),dtype=np.complex) # create the matrix
    if np.max(np.abs(mout.todense() - mout.todense().H))>0.00001: raise
    return mout
  except:
    print "empty matrix"
    mout = csc_matrix((d,d),dtype=np.complex) # create the matrix
    return mout 


def read_orbitals():
  lines = open("orbital.in","r").readlines()
  orb = []
  for l in lines:
    orb.append(l.split()[0]) # append first element
  return orb



def read_basis(namefile):
  m = np.genfromtxt(namefile) # read matrix from file
  from basis import BaseMB
  from copy import deepcopy
  bs = [] # empty list
  orb = read_orbitals() 
  for v in m:
    b = BaseMB() # create object
    b.num2occ(v) # setup the object
    b.orb = orb # orbital atribute
    bs.append(deepcopy(b)) # add to the list
  return bs  # return the list of bassi objects


