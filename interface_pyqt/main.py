import os
import sys
import glob # to find the generated files
import shutil # to copy files in any operative system
import tempfile # to create the temporal folder in any operative system
import subprocess # to call pdflatex and the pdf viewer
mainpath = os.path.dirname(os.path.realpath(__file__))
original_path = os.getcwd() # original path where tranci is being executed
tranciroot = mainpath +"/../" # root to tranci
sys.path.append(tranciroot+"/src/") # add tranci library
sys.path.append(mainpath) # add this folder
import numpy as np
from tranci.read import read_matrix
from tranci.atom import CIatom
from tranci import hamiltonians
from tranci.hamiltonians import eigenvalues,lowest_states
import os # for calling the terminal
from tranci import write # for writing in latex 
import matplotlib.pyplot as py
from tranci.check import check_all # check that the hamiltonian is right
do_check = True # perform check of the hamiltonian

print("Tranci has been executed in",original_path)

from tranci import numberformat
numberformat.tol = 1e-3

import matplotlib.cm as cmplt
import matplotlib
matplotlib.rcParams.update({'font.size': 22}) # increase font size

import qtwrap # import the library with simple wrappers to qt4
get = qtwrap.get  # get the value of a certain variable
getbox = qtwrap.getbox  # get the value of a certain variable
window = qtwrap.main() # this is the main interface
# reject malformed numbers as they are typed, instead of silently zeroing them
qtwrap.add_numeric_validators(["U","n","soc","D","E","trigonal","O","z4",
  "x2y2","B","theta_b","phi_b","j","theta_j","phi_j","tol_ene","nwf_heff",
  "zaxis_x","zaxis_y","zaxis_z","initial_value","final_value","steps",
  "lineEdit"])

## temporal folder, in the location for temporal files of this system
temporal_folder = tempfile.mkdtemp(prefix="tranci_tmp_") # portable temporal folder

# go to a temporal folder
def restart_tranci():
  """Fully restart tranci"""
  tmpfol = temporal_folder
  print("Temporal Tranci folder is",tmpfol)
  os.chdir(tmpfol) # go to the temporal folder

def save_tranci():
    """Save the generated data in a file"""
    save_folder = os.path.join(original_path,"tranci_data") # name of the folder
    def tcopy(a):
        for f in glob.glob(os.path.join(temporal_folder,a)): # loop over matches
            shutil.copy(f,save_folder) # copy this file
    os.makedirs(save_folder,exist_ok=True) # create the folder if not there
    tcopy("*.tex")
    tcopy("*.pdf")
    tcopy("*.OUT")
    print("Saved tranci data in ",save_folder)



restart_tranci() # restart tranci


def show_pdf():
    name = os.path.abspath("spectrum_ci.pdf") # full path to the pdf
    if not os.path.isfile(name): # no pdf generated yet
      print("No pdf found, run the calculation first")
      return
    if sys.platform=="win32": os.startfile(name) # Windows system
    elif sys.platform=="darwin": subprocess.Popen(["open",name]) # Mac system
    else: subprocess.Popen(["xdg-open",name]) # Linux system



def run_pdflatex():
  """Compile the latex summary, returning True on success"""
  r = subprocess.run(["pdflatex","-interaction=nonstopmode","spectrum_ci.tex"],
          stdout=subprocess.PIPE,stderr=subprocess.STDOUT) # capture the output
  if r.returncode!=0: # show why, instead of claiming success later
    out = r.stdout.decode("utf-8","replace") if r.stdout else ""
    print("pdflatex failed (exit %d). Last lines of its output:"%r.returncode)
    for l in out.splitlines()[-20:]: print("  "+l)
    return False
  return True



def sweep_label(stype):
  """Axis label for a swept variable

  The angular entries of the sweep combobox are theta_B/phi_B/theta_J/phi_J,
  and get_b multiplies them by pi, so they are in units of pi radians."""
  if stype in ["theta_B","phi_B","theta_J","phi_J"]: return stype+" [$\\pi$ rad]"
  return stype+" [eV]"


def get_b(babs,theta,phi):
  """Get the magnetic field"""
  st = np.sin(theta*np.pi)
  ct = np.cos(theta*np.pi)
  sp = np.sin(phi*np.pi)
  cp = np.cos(phi*np.pi)
  b = babs*np.array([st*cp,st*sp,ct])    # build the magnetic field
  return b



def read_inputs():
  """Read all the inputs"""
  class params: pass
  p = params() # parameters of the system
  p.n = int(get("n"))
  p.D = get("D")
  p.z4 = get("z4")
  p.E = get("E")
  p.x2y2 = get("x2y2")
  p.O = get("O")
  p.U = get("U")
  p.trigonal = get("trigonal")
  p.soc = get("soc")
  p.theta_b = get("theta_b") # theta for magnetic field
  p.phi_b = get("phi_b") # phi for Zeeman
  p.theta_j = get("theta_j") # theta for magnetic field
  p.phi_j = get("phi_j") # phi for Zeeman
  p.babs = get("B") # absolute Zeeman
  p.jabs = get("j") # absolute Zeeman
  p.b = get_b(p.babs,p.theta_b,p.phi_b) # get the magnetic field
  p.j = get_b(p.jabs,p.theta_j,p.phi_j) # get the magnetic field
  return p

def initialize_one_shot():
  """ Initialize the one shot calculation"""
  p = read_inputs() # read all the inputs
#  os.system("cp "+ str(p.n) + "/* ./") # copy input files
  at = get_atom() # read the basis from file
  at.wavefunction_z_axis = [get("zaxis_x"),get("zaxis_y"),get("zaxis_z")]
  m = hamiltonians.build_hamiltonian(at,p) # get the hamiltonian
  if do_check:  check_all(at) # check the hamiltonian
  header = hamiltonians.latex_DE(at,p) # string for the hamiltonian
  ls = lowest_states(m,atom=at) # create the object
  ls.disentangle_manifolds(at.jz) # disentangle manifold
#  ls.get_gtensor() # compute gfactor
#  exit()
  if get("nwf_heff")>1: nw = int(get("nwf_heff"))
  else: nw=None
  write.write_all(ls,header=header,n=nw) # write Latex file
  if shutil.which("pdflatex") is None: # no latex in this system
    print("pdflatex not found, only spectrum_ci.tex was written")
    print("Install a LaTeX distribution (MiKTeX in Windows) for the pdf summary")
    return
  ok = run_pdflatex() # compile the latex file
  ok = run_pdflatex() and ok # do it twice, so that the index is right
  if not ok or not os.path.isfile("spectrum_ci.pdf"): # do not claim success
    print("The PDF summary could not be created; spectrum_ci.tex was written")
    return
  print("#########################")
  print("## PDF Summary created ##")
  print("#########################")



def get_atom():
  p = read_inputs() # read all the inputs
#  os.system("cp "+ str(p.n) + "/* ./") # copy input files
#  os.system("cp "+tranciroot+"cilib/"+ str(p.n) + "/* ./") # copy input files
  from tranci import atom
  at = atom.get_atom(ne=p.n)
#  at = CIatom() # create the CI object
#  at.read() # read all the matrices
#  at.get_basis() # read the basis from file
  hamiltonians.tol = np.max([1e-8,get("tol_ene")]) 
  hamiltonians.ntol = -int(round(np.log10(hamiltonians.tol)))
  return at # return atom




def initialize_sweep():
  """Launch a sweeping calculation"""
  p = read_inputs() # read all the inputs
#  os.system("cp "+ str(p.n) + "/* ./") # copy input files
  at = get_atom()
  def fsweep(x):
    """Function to perform the sweep"""
    stype = getbox("sweep_variable") # get the variable
    if stype == "U": p.U = x
    elif stype == "z^2": p.D = x
    elif stype == "x^2-y^2": p.E = x
    elif stype == "soc": p.soc = x
    elif stype == "z^4": p.z4 = x
    elif stype == "x^2y^2": p.x2y2 = x
    elif stype == "x^4+y^4+z^4": p.O = x
    elif stype == "(x+y+z)^2": p.trigonal = x
    elif stype == "B": p.b = get_b(x,p.theta_b,p.phi_b)
    elif stype == "theta_B": p.b = get_b(p.babs,x,p.phi_b)
    elif stype == "phi_B": p.b = get_b(p.babs,p.theta_b,x)
    elif stype == "J": p.j = get_b(x,p.theta_j,p.phi_j)
    elif stype == "theta_J": p.j = get_b(p.jabs,x,p.phi_j)
    elif stype == "phi_J": p.j = get_b(p.jabs,p.theta_j,x)
    else: raise # raise error
    m = hamiltonians.build_hamiltonian(at,p) # get the hamiltonian
    ls = lowest_states(m,atom=at) # perform the calculation 
    return ls # return the object
  return fsweep # return function

def plot_eigenvalues(write=True,center=True):
  """Plots the excited states"""
  ###############################
  ###############################
  fsweep = initialize_sweep() # get the generator function
  xs = get_sweep_parameters() # get the array
  gst = [fsweep(ix) for ix in xs]  # create the list of objects
#  ys = [g.get_excitations() for g in gst] # get the energies of the excitations
  ys = [g.evals_full for g in gst] # get all the eigenvalues
  fig = py.figure() # create figure
  fig.subplots_adjust(.2,.15) # adjust the subplots
  ys = np.array(ys).transpose() # row is same eigenvector evolving
  # number of energies to plot
  nenergies = int(get("lineEdit")) # number of energies to plot
  if 0 < nenergies < len(ys): ys = np.array([ys[i] for i in range(nenergies)])
  else: pass

  # now move the center of gravity to wherever it should be
  if center: 
    ys = ys.transpose() # row is same value
    ys = [y - sum(y)/len(y) for y in ys] # with respect to the center
    ys = np.array(ys).transpose() # transpose, row is eigenvalue evolving
  else: 
    ys = [y - ys[0] for y in ys] # with respect GS
    ys = np.array(ys) # convert to array, row is eigenvalue evolving
  # now plot
  colors = cmplt.rainbow(np.linspace(0, 1, len(ys))) # different colors
  fo = open("EIGENVALUES.OUT","w")
  for (y,c) in zip(ys,colors): # loop over eigenvalue
    py.plot(xs,y,marker="o",c=c) 
    for (ix,iy) in zip(xs,y): fo.write(str(ix)+" "+str(iy)+"\n")
  fo.close()
  py.xlim([min(xs),max(xs)]) # 
  py.ylabel("Energy [eV]")  # label for the y axis
  stype = getbox("sweep_variable")
  py.xlabel(sweep_label(stype))  # label for the x axis
  fig.set_facecolor("white")
  py.tight_layout()
  py.show()




def plot_spectrum():
  """ Plot the different eigenvlaues, centered"""
  plot_eigenvalues(center=True)


def plot_excitations():
  """ Plot the different eigenvlues, shifted to the GS"""
  plot_eigenvalues(center=False)

def get_sweep_parameters():
  return np.linspace(get("initial_value"),get("final_value"),int(get("steps"))) 


def plot_degeneracy():
  """Plots the degeneracy of the ground state"""
  ###############################
  ###############################
  xs = get_sweep_parameters() # get the array
  fsweep = initialize_sweep() # get the generator function
  gst = [fsweep(ix) for ix in xs]  # create the list of objects
  T = hamiltonians.tol # tolerancy
  ds = [g.get_gs_multiplicity(tol=T) for g in gst] # get degeneracies
  fig = py.figure() # create figure
  fig.subplots_adjust(.2,.15) # adjust the subplots
  py.plot(xs,ds,c="green",marker="o") 
  np.savetxt("DEGENERACY.OUT",np.array([xs,ds]).T) # save data
  py.xlim([min(xs),max(xs)]) # 
  py.ylim([0,max(ds)+1]) # 
  py.ylabel("Degeneracy")  # label for the y axis
  stype = getbox("sweep_variable")
  py.xlabel(sweep_label(stype))  # label for the x axis
  fig.set_facecolor("white")
  py.tight_layout()
  py.show()





def plot_operator():
  """Plots the eigenvalues of a certain operator in the GS"""
  ###############################
  ###############################
  fsweep = initialize_sweep() # get the generator function
  xs = get_sweep_parameters() # get the array
  gst = [fsweep(ix) for ix in xs]  # create the list of objects
  at = get_atom() # get the atom object

  ########################
  ########################
  ########################
  oname = getbox("operator_name")
  if (oname=="Sx"): op = at.sx     # get this operator 
  elif (oname=="Sy"): op = at.sy   # get this operator 
  elif (oname=="Sz"): op = at.sz   # get this operator 
  elif (oname=="S_(111)"): op = (at.sz + at.sx + at.sy)/np.sqrt(3)
  elif (oname=="L_(111)"): op = (at.lz + at.lx + at.ly)/np.sqrt(3)
  elif (oname=="Jx"): op = at.jx   # get this operator
  elif (oname=="Jy"): op = at.jy   # get this operator
  elif (oname=="Jz"): op = at.jz   # get this operator
  elif (oname=="Lx"): op = at.lx   # get this operator 
  elif (oname=="Ly"): op = at.ly   # get this operator 
  elif (oname=="Lz"): op = at.lz   # get this operator 
  elif (oname=="L2"): op = at.l2   # get this operator 
  elif (oname=="S2"): op = at.s2   # get this operator 
  elif (oname=="J2"): op = at.j2   # get this operator 
  elif (oname=="x2"): op = at.x2   # get this operator 
  elif (oname=="y2"): op = at.y2   # get this operator 
  elif (oname=="z2"): op = at.z2   # get this operator 
  elif (oname=="LS"): op = at.ls   # get this operator 
  elif (oname=="up m=-2"): op = at.um2   # get this operator 
  elif (oname=="up m=-1"): op = at.um1   # get this operator 
  elif (oname=="up m=0"): op = at.u0   # get this operator 
  elif (oname=="up m=+1"): op = at.up1   # get this operator 
  elif (oname=="up m=+2"): op = at.up2   # get this operator 
  elif (oname=="dn m=+2"): op = at.dm2   # get this operator 
  elif (oname=="dn m=-1"): op = at.dm1   # get this operator 
  elif (oname=="dn m=0"): op = at.d0   # get this operator 
  elif (oname=="dn m=+1"): op = at.dp1   # get this operator 
  elif (oname=="dn m=+2"): op = at.dp2   # get this operator 
  elif (oname=="x4y4z4"): op = at.x4+at.y4+at.z4   # get this operator 
  else: print(oname) ; raise
  ########################
  ########################
  ########################

  evals = [g.get_gs_projected_eigenvalues(op) for g in gst] # get op eigen 
  fig = py.figure() # create figure
  fig.subplots_adjust(.2,.15) # adjust the subplots
  #############
  fo = open("OPERATOR_VALUES.OUT","w")
  for (x,y) in zip(xs,evals):
    colors = cmplt.rainbow(np.linspace(0, 1, len(y))) # different colors
    xplot = [x for iy in y]
    py.scatter(xplot,y,c=colors) 
    for (ix,iy) in zip(xplot,y): fo.write(str(ix)+" "+str(iy)+"\n")
  fo.close() # close file
  py.xlim([min(xs),max(xs)]) # 
#  py.ylim([min(evals),max(evals)]) # 
  py.ylabel("$\\langle "+oname+"\\rangle$")  # label for the y axis
  stype = getbox("sweep_variable")
  py.xlabel(sweep_label(stype))  # label for the x axis
  fig.set_facecolor("white")
  py.tight_layout()
  py.show() # show graph













# Add the tranci logo
from qtwrap import set_logo
tranci_logo = tranciroot+"/logos/orbitals.png"
set_logo("logo",tranci_logo)


# create signals
signals = dict()
#signals["on_window_destroy"] = gtk.main_quit  # close the window
signals["initialize_one_shot"] = initialize_one_shot  # initialize and run
signals["plot_spectrum"] = plot_spectrum  # initialize and run
signals["plot_excitations"] = plot_excitations  # initialize and run
signals["plot_degeneracy"] = plot_degeneracy  # initialize and run
signals["plot_operator"] = plot_operator  # initialize and run
signals["show_pdf"] = show_pdf  # show pdf with the results
signals["save_tranci"] = save_tranci  # show pdf with the results

window.connect_clicks(signals) 
window.run()

