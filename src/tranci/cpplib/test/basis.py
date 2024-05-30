import numpy as np

class BaseMB():
  """Class for a many body base vector"""
  def num2occ(self,v):
    """ Sets up the occupation of the vector using as input an array of 0,1"""
    self.occ = [] # empty list
    for o in v: # 
      if int(o)==1:
        self.occ.append(True)
      elif int(o)==0:
        self.occ.append(False)
      else:
        raise
  def get_latex(self):
    """ Return the current vector in latex form """
#    l = "" # empty string
    l = "|" # empty string
    first = True
    for (i,orb) in zip(range(len(self.occ)),self.occ):
      if orb:
#        l += "c^\dagger_{"+self.orb[i]+"}"
        if not first: l+=" , "
        else:  first = False          
        l += self.orb[i]
#    l += "|\Omega\\rangle"
    l += "\\rangle"
    return l # return the latex line







