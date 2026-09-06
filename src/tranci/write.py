# functions that write the LaTeX summary of a calculation

import time
import numpy as np
import scipy.linalg as lg
from .numberformat import zform
from .numberformat import texnum
from .numberformat import recognise_number
from . import templates
from . import hamiltonians


# how many excited manifolds are written in the body of the document, the
# rest go to an appendix so that the main text stays short
nmanifolds_body = 10


PREAMBLE = r"""\documentclass[11pt,a4paper]{article}
\usepackage[T1]{fontenc}
\usepackage[utf8]{inputenc}
\usepackage{lmodern}
\usepackage{microtype}
\usepackage{amsmath}
\usepackage{amssymb}
\usepackage{array}
\usepackage{booktabs}
\usepackage{longtable}
\usepackage{xcolor}
\usepackage[margin=2.5cm]{geometry}
\usepackage{parskip}
\usepackage[colorlinks=true,linkcolor=blue!50!black,urlcolor=blue!50!black]{hyperref}
\setcounter{tocdepth}{1}
\setlength{\emergencystretch}{3em}
\allowdisplaybreaks
\newcommand{\ket}[1]{\left|#1\right\rangle}
"""


def build_latex(formula,name,title="Configuration interaction summary"):
  """Build a self-contained LaTeX document around a body string"""
  stamp = time.strftime("%Y-%m-%d %H:%M")
  f = open(name+".tex","w")
  f.write(PREAMBLE)
  f.write("\\title{\\textbf{TranCI}\\\\[0.4em]\\large "+title+"}\n")
  f.write("\\author{}\n")
  f.write("\\date{"+stamp+"}\n")
  f.write("\\begin{document}\n")
  f.write("\\maketitle\n")
  f.write("\\tableofcontents\n\\newpage\n\n")
  f.write(formula)
  f.write("\n\\end{document}\n")
  f.close() # close the file


def longtable(align,header,rows,caption=None):
  """Return a booktabs longtable that breaks across pages on its own.
  align is a column specification, header a list of column titles and
  rows a list of lists of already formatted cells"""
  nc = len(header)
  head = "  "+" & ".join(header)+" \\\\\n"
  out = "\\begin{center}\n\\begin{longtable}{"+align+"}\n"
  if caption is not None: out += "\\caption{"+caption+"}\\\\\n"
  out += "\\toprule\n"+head+"\\midrule\n\\endfirsthead\n"
  out += "\\toprule\n"+head+"\\midrule\n\\endhead\n"
  out += "\\midrule\n\\multicolumn{"+str(nc)+"}{r}{\\small\\emph{continued on next page}}\\\\\n\\endfoot\n"
  out += "\\bottomrule\n\\endlastfoot\n"
  for r in rows: out += "  "+" & ".join(r)+" \\\\\n"
  out += "\\end{longtable}\n\\end{center}\n\n"
  return out


def get_energies_table(lowest):
  """Return the spectrum section: a table of manifolds, degeneracies and
  energies"""
  dg = lowest.get_degeneracies() # get the degeneracies (# and energy)
  ex = lowest.get_excitations() # get the excitations
  numconf = sum([d[0] for d in dg])
  text = "\\section{Spectrum}\n\n"
  text += "The Hamiltonian has "+str(numconf)+" eigenstates, grouped into "
  text += str(len(dg))+" manifolds of degenerate levels. "
  e0 = getattr(lowest,"e0",None) # absolute ground state energy, if known
  if e0 is not None:
    text += "The ground state energy is $E_0 = "+texnum(e0,6)+"$ eV, and "
    text += "every energy below is measured from it.\n\n"
  else: text += "Energies are measured from the ground state.\n\n"
  header = ["Manifold","Degeneracy","$\\Delta E$ (eV)","$\\Delta E$ (meV)"]
  rows = []
  for (i,d) in enumerate(dg): # loop over manifolds
     name = "Ground state" if i==0 else "Excited "+str(i)
     rows.append([name,str(d[0]),texnum(ex[i],6),texnum(ex[i]*1000,3)])
  text += longtable("lcrr",header,rows)
  return text


def get_table_states(lowest,ghz=True):
  """Return the section with the energy and expectation values of every
  eigenstate"""
  at = lowest.atom # get the object
  text = "\\section{Eigenstates and expectation values}\n\n"
  text += "Energy and expectation values of the orbital and spin angular "
  text += "momentum of every eigenstate, ordered by increasing energy. "
  text += "Energies are measured from the ground state.\n\n"
  m = templates.write_low_energy(lowest,at,n=len(lowest.evals),ghz=ghz)
  m = np.atleast_2d(m) # a single state comes back as a 1d array
  e0 = m[0][0],m[0][1]
  header = ["$\\Delta E$ (GHz)","$\\Delta E$ (meV)","$\\langle L_x\\rangle$",
            "$\\langle L_y\\rangle$","$\\langle L_z\\rangle$",
            "$\\langle S_x\\rangle$","$\\langle S_y\\rangle$",
            "$\\langle S_z\\rangle$"]
  rows = []
  for im in m: # loop over states
     row = [texnum(im[0]-e0[0],3),texnum(im[1]-e0[1],3)]
     row += [texnum(x,3) for x in im[2:]]
     rows.append(row)
  text += longtable("rrrrrrrr",header,rows)
  return text


def write_energies(lowest):
  """ Write a table with the energies and degeneracies"""
  text = get_energies_table(lowest)
  build_latex(text,"formula")


def get_gs_operators(lowest):
  """Return the section with the operators projected onto the ground state
  manifold"""
  def gr(m):
    """ Return the representation of a certain matrix"""
    return hamiltonians.get_representation(lowest.gs_manifold,m)
  atom = lowest.atom # get object
  ops = [(atom.sx,"S_x"),(atom.sy,"S_y"),(atom.sz,"S_z"),(atom.s2,"S^2"),
         (atom.lx,"L_x"),(atom.ly,"L_y"),(atom.lz,"L_z"),(atom.l2,"L^2"),
         (atom.jx,"J_x"),(atom.jy,"J_y"),(atom.jz,"J_z"),(atom.j2,"J^2"),
         (atom.x2,"\\sum_i l_{x,i}^2"),(atom.y2,"\\sum_i l_{y,i}^2"),
         (atom.z2,"\\sum_i l_{z,i}^2"),
         (atom.x4+atom.y4+atom.z4,"\\sum_i (l_{x,i}^4+l_{y,i}^4+l_{z,i}^4)"),
         (atom.vc,"V_{e\\text{-}e}"),
         (atom.ls,"\\sum_i \\vec l_i \\cdot \\vec s_i")]
  ndeg = len(lowest.gs_manifold)
  text = "\\section{Operators projected onto the ground state}\n\n"
  if ndeg==1: # a single state: expectation values only
    text += "The ground state is non-degenerate, so every operator reduces "
    text += "to its expectation value $\\langle\\Psi_1|O|\\Psi_1\\rangle$.\n\n"
    rows = []
    for (m,name) in ops:
      v = gr(m)[0,0]
      s = texnum(v,4)
      if abs(np.imag(v))>1e-6: s += " $"+zform(v)+"$"
      rows.append(["$"+name+"$",s])
    text += longtable("lr",["Operator","Expectation value"],rows)
    return text
  text += "The ground state manifold has "+str(ndeg)+" states, and every "
  text += "operator is written as a matrix in that manifold. The eigenvalues "
  text += "of each block are listed below it.\n\n"
  for (m,name) in ops: text += matrix2latex(gr(m),name)
  return text


def write_gs_operators(atom,lowest):
  """Writes several operators in a file"""
  text = get_gs_operators(lowest)
  build_latex(text,"formula")


def matrix2latex(m,name=""):
  """Return a matrix as a displayed equation, followed by its eigenvalues"""
  m = np.array(m)
  n = len(m) # length of the matrix
  lhs = name+" = " if name!="" else ""
  if n>10:
    return "\n"+lhs.replace(" = ","")+" is a "+str(n)+"$\\times$"+str(n)+\
           " matrix, too large to be displayed.\n\n"
  mstr = "\\begin{equation*}\n"+lhs+"\\begin{pmatrix}\n"
  for i in range(n):
    mstr += "  "+" & ".join([zform(m[i,j]) for j in range(n)])+" \\\\\n"
  mstr += "\\end{pmatrix}\n\\end{equation*}\n"
  es = lg.eigvalsh(m) # compute eigenvalues
  mstr += "\\noindent Eigenvalues: $\\{"
  mstr += ",\\; ".join([zform(e) for e in es])
  mstr += "\\}$\n\n"
  return mstr


def format_wavefunctions(waves):
  """Format a list of wavefunctions, each a list of (coefficient, ket)
  strings, as displayed equations"""
  form = ""
  for (ii,v) in enumerate(waves): # loop over wavefunctions
    form += "\\begin{equation*}\n\\Psi_{"+str(ii+1)+"} = "
    form += "\\left(\\begin{array}{r@{\\;\\;}l}\n"
    for iv in v: form += "  "+iv+" \\\\\n"
    form += "\\end{array}\\right)\n\\end{equation*}\n\n"
  return form


def manifold_block(at,lowest,i):
  """Return the heading and wavefunctions of the i-th manifold"""
  dg = lowest.get_degeneracies()
  ex = lowest.get_excitations()
  vs = lowest.manifolds[i]
  vs = [at.rotate_wavefunction_axis(v) for v in vs] # rotate wavefunctions
  formulas = [at.get_latex_wavefunction(v) for v in vs]
  if i==0: title = "Ground state"
  else: title = "Excited manifold "+str(i)
  title += " ($\\Delta E = "+texnum(ex[i]*1000,3)+"$ meV, "
  title += str(len(vs))+(" state" if len(vs)==1 else " states")+")"
  return "\\subsection*{"+title+"}\n\n"+format_wavefunctions(formulas)


def get_manifolds(at,lowest,nbody=None):
  """Return the wavefunction sections. The ground state and the first nbody
  excited manifolds go in the main text, the rest to an appendix"""
  if nbody is None: nbody = nmanifolds_body
  z = at.wavefunction_z_axis
  text = "\\section{Wavefunctions}\n\n"
  text += "Eigenstates written in the basis of Slater determinants "
  text += "$\\ket{m_1\\sigma_1,\\,m_2\\sigma_2,\\ldots}$, with $m$ the orbital "
  text += "angular momentum projection and $\\sigma$ the spin along the "
  text += "quantization axis $\\hat z = ("
  text += ",\\,".join([recognise_number(zi) for zi in z])+")$. "
  text += "Only components with weight above 1\\% are shown, sorted by "
  text += "decreasing weight; a missing coefficient means the state is a "
  text += "single determinant.\n\n"
  nm = len(lowest.manifolds)
  nshow = min(nm,nbody+1)
  for i in range(nshow): text += manifold_block(at,lowest,i)
  if nshow<nm:
    text += "The remaining "+str(nm-nshow)+" manifolds are listed in "
    text += "Appendix~\\ref{app:manifolds}.\n\n"
  return text


def get_manifolds_appendix(at,lowest,nbody=None):
  """Return the appendix with the manifolds not shown in the main text"""
  if nbody is None: nbody = nmanifolds_body
  nm = len(lowest.manifolds)
  nshow = min(nm,nbody+1)
  if nshow>=nm: return ""
  text = "\\appendix\n\\section{Higher excited manifolds}\n"
  text += "\\label{app:manifolds}\n\n"
  for i in range(nshow,nm): text += manifold_block(at,lowest,i)
  return text


def write_manifolds(at,lowest):
  """ Writes the manifolds in a file"""
  text = get_manifolds(at,lowest)+get_manifolds_appendix(at,lowest)
  build_latex(text,"formula")


def write_all(lowest,header="",n=None,nbody=None):
  """Write the full LaTeX summary to spectrum_ci.tex. header is the
  Hamiltonian section, n the number of states fitted with an effective
  Hamiltonian (None to skip it) and nbody the number of excited manifolds
  written in the main text before the appendix"""
  at = lowest.atom
  # these populate lowest.gs_manifold / lowest.manifolds, which the helpers
  # below read directly; without them a library caller got an AttributeError
  if not hasattr(lowest,"gs_manifold"): lowest.get_gs_manifold()
  if not hasattr(lowest,"manifolds"): lowest.get_manifolds()
  text = header
  text += get_energies_table(lowest)
  text += get_table_states(lowest)
  text += get_gs_operators(lowest)
  text += get_manifolds(at,lowest,nbody=nbody)
  text += get_effective_hamiltonian(lowest,n=n)
  text += get_manifolds_appendix(at,lowest,nbody=nbody)
  build_latex(text,"spectrum_ci")  # name of the file


def get_effective_hamiltonian(lowest,n=None):
    """Return latex form of the effective Hamiltonian"""
    if n is None: return ""
    try: from .effectivehamiltonian import effective_hamiltonian
    except ImportError: # jax is not installed in this system
        print("jax not found, skipping the effective Hamiltonian")
        return "\\section{Effective Hamiltonian}\n\nOmitted: jax is not installed.\n\n"
    text = effective_hamiltonian(lowest,n=n) # return the effective Hamiltonian
    return text + "\n\n"


def write_gs_manifold(at,lowest):
  """Writes the ground state manifold in a file"""
  vs = lowest.gs_manifold # ground state manifold
  formulas = [at.get_latex_wavefunction(v) for v in vs]
  form = format_wavefunctions(formulas)
  build_latex(form,"formula")
