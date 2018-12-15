import numpy as np
from read import read_matrix
from atom import CIatom
from hamiltonians import DE,eigenvalues,lowest_states,latex_DE

##########################
##  Begining of inputs
##########################
D = 0.02 # uniaxial crystal field
E = 0.0 # 
O = 0.0 # octahedral crystal field
U = 2.0 # electron electron interaction
n = 2 # number of electrons
soc = 0. # spin orbit coupling
z4 = 0.02 # cuartic axial symmetry
babs = 0.0
theta = 0.0
phi = 0.0

##########################
##  End of inputs
##########################


import os
for n in range(1,10):
  os.system("cp ../data/"+ str(n) + "/* ./") # copy input files
  
  #################################
  #################################
  #################################
  at = CIatom() # create the CI object
  at.read() # read all the matrices
  at.get_basis() # read the basis from file
  
  import check
  try:
    check.check_all(at)
    print "Passed ",n
  except:
    print "Error in",n
    break  

