# this functions write different latex documents

from format import fform
from format import zform

import hamiltonians

def get_gs_operators(atom,lowest):
  """Gets the latex representation of several ground state operators"""
  def gr(m): 
    """ Return the representation of a certain matrix"""
    return hamiltonians.get_representation(lowest.gs_manifold,m) # get Lz
  mstr = "\\section{Proyection of operator on GS}\n"
  mstr += matrix2latex(gr(atom.sx),"S_x") # get the matrix in latex
  mstr += matrix2latex(gr(atom.sy),"S_y") # get the matrix in latex
  mstr += matrix2latex(gr(atom.sz),"S_z") # get the matrix in latex
  mstr += matrix2latex(gr(atom.s2),"S^2") # get the matrix in latex
  mstr += matrix2latex(gr(atom.lx),"L_x") # get the matrix in latex
  mstr += matrix2latex(gr(atom.ly),"L_y") # get the matrix in latex
  mstr += matrix2latex(gr(atom.lz),"L_z") # get the matrix in latex
  mstr += matrix2latex(gr(atom.l2),"L^2") # get the matrix in latex
  mstr += matrix2latex(gr(atom.jx),"J_x") # get the matrix in latex
  mstr += matrix2latex(gr(atom.jy),"J_y") # get the matrix in latex
  mstr += matrix2latex(gr(atom.jz),"J_z") # get the matrix in latex
  mstr += matrix2latex(gr(atom.j2),"J^2") # get the matrix in latex
  mstr += matrix2latex(gr(atom.x2),"X^2") # get the matrix in latex
  mstr += matrix2latex(gr(atom.y2),"Y^2") # get the matrix in latex
  mstr += matrix2latex(gr(atom.z2),"Z^2") # get the matrix in latex
  mstr += matrix2latex(gr(atom.x4+atom.y4+atom.z4),"X^4+Y^4+Z^4") # get the matrix in latex
  mstr += matrix2latex(gr(atom.vc),"V_{e-e}") # get the matrix in latex
  mstr += matrix2latex(gr(atom.ls),"\\vec L^{(i)} \cdot \\vec S^{(i)}") # get the matrix in latex
  return mstr

def write_gs_operators(atom,lowest):
  """Writes several operators in a file"""
  text = get_gs_operators(atom,lowest)
  build_latex(text,"formula")


def get_energies_table(lowest):
  """ Returns a string with a latex form of the degeneracies"""
  table = "" # initialice the table
  table += "\\section{Spectrum}\n" # initialice the table
  beginning = "" # string for beggining the table
  beginning += "\\begin{center}\n"
  beginning += "\\begin{tabular}{| c | c | c| c |}\n"
  beginning += "\\hline\n" # horizontal line
  beginning +="Manifold  &  Degeneracy  & $\Delta E$ (eV) & Energy (eV) "+"\\"+"\\"+"\n"
  beginning += "\\hline\n" # horizontal line
  end = "" # string for ending the table
  end += "\\hline\n" # horizontal line
  end += "\\end{tabular}\n"
  end += "\\end{center}\n"
  dg = lowest.get_degeneracies() # get the degeneracies (# and energy)
  ex = lowest.get_excitations() # get the excitations 
  istate = 0 # counter
  table += beginning  # beginning of the table
  for d in dg: # loop over manifolds
     if istate==0: name = "GS" # name for ground state
     else: name = "Excited \#"+str(istate) # name for excited
     table += name + "  &  "  # line with name 
     table += str(d[0]) +"   &  " # line with degeneracy
     table += fform(ex[istate]) +"  &  " # line with distance to ground
     table += fform(d[1]) # line with energy 
     table += "  \\"+"\\"+"\n"  # line with name degeneracy and energy
     istate += 1
     if istate%36 ==0: table += end + "\n\n" + beginning # new page
  table += end # end of the table
  return table

def write_energies(lowest):
  """ Write a table with the energies and degeneracies"""
  text = get_energies_table(lowest)
  build_latex(text,"formula")


def get_manifolds(at,lowest):
  """Writes the ground state manifold in a file"""
  allform = "" # empty string for the total latex stuff
  allform += "\\section{Wavefunctions}\n"
  istate = 0 # counter for the manifold
  for vs in lowest.manifolds: # loop over manifolds
    formulas = [] # list for the formulas
    for v in vs:
      formulas.append(at.get_latex_wavefunction(v)) # get the wavefunction
    form = format_wavefunctions(formulas)
    if istate==0:   title = "\n\n\subsection{Ground state manifold}\n\n"
    else:   title = "\n\n\subsection{Excited state manifold \#"+str(istate)+"}\n\n"
    allform += title + form # append this string
    istate += 1
  return allform # add all the formulas


def write_manifolds(at,lowest):
  """ Writes the manifolds in a file"""
  text = get_manifolds(at,lowest)
  build_latex(text,"formula")



def write_all(at,lowest,header=""):
  """ Writes all the stuff in a file"""
  text = header
  text += get_energies_table(lowest)
  text += get_gs_operators(at,lowest)
  text += get_manifolds(at,lowest)
  build_latex(text,"spectrum_ci")  # name of the file



def matrix2latex(m,name=""):
  """ Returns a matrix in latex form """
  mstr = "\n" # initialize matrix string
  n = len(m) # length of the matrix
  if n>20: 
    mstr += "\\text{Matrix too large..}"
    return mstr # skip
  mstr += "\\begin{equation}\n"
  mstr += name + " = \n" # name of the matrix
  mstr += "\\begin{pmatrix}\n"
  for i in range(n):
    for j in range(n):
      mstr += zform(m[i,j])
      if j<(n-1): mstr += "  &  " # not last number
      else: mstr += "  \\"+"\\"+"\n" # last number       
  mstr += "\\end{pmatrix}\n"
  mstr += "\\end{equation}\n"
  return mstr





def write_gs_manifold(at,lowest):
  """Writes the ground state manifold in a file"""
  vs = lowest.gs_manifold # ground state manifold
  formulas = [] # list for the formulas
  for v in vs:
#    print at.get_latex_wavefunction(v) # get the wavefunction
    formulas.append(at.get_latex_wavefunction(v)) # get the wavefunction
  form = format_wavefunctions(formulas)
  build_latex(form,"formula")


def write_gs_projected_operator(at,ls,operators):
  """Writes the projection of certain operators on the ground state manifold"""
  vs = ls.gs_manifold # ground state manifold
  










def format_wavefunctions(waves):
  """Formats a list of wavefunctions to be written in latex"""
  ii = 0 # counter
  form="\n"
  for v in waves: # loop over wavefunctions
    form += "\\begin{equation}\n" # output string
    form += "\Psi_{"+str(ii+1) + "} = \n" # add 
    ii += 1 # increase counter
    form += "\\begin{pmatrix}\n" # output string
    for iv in v:
      form += iv +"  \\"+"\\"+"\n"
    form += "\\end{pmatrix}\n" # output string
    form += "\\end{equation}\n\n" # output string
  return form    


def build_latex(formula,name):
  """ Build a latex page with a certain formula"""
  f = open(name+".tex","w")
  f.write("\documentclass{article}\n")
  f.write("\usepackage{amsmath}\n")
  f.write("\\begin{document}\n")
  f.write(formula)
  f.write("\\end{document}\n")
  f.close() # close the file




