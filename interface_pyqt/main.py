import os
import sys
import glob # to find the generated files
import json # to save and load the parameters
import shutil # to copy files in any operative system
import tempfile # to create the temporal folder in any operative system
import subprocess # to call pdflatex and the pdf viewer
from contextlib import contextmanager
mainpath = os.path.dirname(os.path.realpath(__file__))
original_path = os.getcwd() # original path where tranci is being executed
tranciroot = mainpath +"/../" # root to tranci
sys.path.append(tranciroot+"/src/") # add tranci library
sys.path.append(mainpath) # add this folder
import numpy as np
from tranci import hamiltonians
from tranci.hamiltonians import lowest_states
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

from PyQt5 import QtCore,QtWidgets
import qtwrap # import the library with simple wrappers to PyQt5
get = qtwrap.get  # get the value of a certain variable
getbox = qtwrap.getbox  # get the value of a certain variable
status = qtwrap.status # message in the status bar
window = qtwrap.main() # this is the main interface
# reject malformed numbers as they are typed, instead of silently zeroing them
# (integer fields are QSpinBox widgets and validate themselves)
qtwrap.add_numeric_validators(["U","soc","D","E","trigonal","O","z4",
  "x2y2","B","theta_b","phi_b","j","theta_j","phi_j","tol_ene",
  "zaxis_x","zaxis_y","zaxis_z","initial_value","final_value"])

# progress bar in the status bar, shown only during sweeps
progress = QtWidgets.QProgressBar()
progress.setMaximumWidth(200)
progress.setTextVisible(True)
window.statusbar.addPermanentWidget(progress)
progress.hide()

parameters_file = "parameters.json" # written next to the results of each run

## temporal folder, in the location for temporal files of this system
temporal_folder = tempfile.mkdtemp(prefix="tranci_tmp_") # portable temporal folder

# go to a temporal folder
def restart_tranci():
  """Fully restart tranci"""
  tmpfol = temporal_folder
  print("Temporal Tranci folder is",tmpfol)
  os.chdir(tmpfol) # go to the temporal folder


@contextmanager
def busy(message):
  """Block the interface, show a wait cursor and a status message while
  a calculation runs; everything is restored even if it raises"""
  status(message)
  QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
  window.centralwidget.setEnabled(False) # no re-entrant clicks
  window.menubar.setEnabled(False)
  qtwrap.app.processEvents()
  try: yield
  except Exception:
    status("Failed: "+message.rstrip(".")) # the error dialog follows
    raise
  finally:
    progress.hide()
    window.centralwidget.setEnabled(True)
    window.menubar.setEnabled(True)
    QtWidgets.QApplication.restoreOverrideCursor()
    qtwrap.app.processEvents()


def save_tranci():
    """Save the generated data in a file"""
    save_folder = os.path.join(original_path,"tranci_data") # name of the folder
    patterns = ["*.tex","*.pdf","*.OUT","*.json"]
    files = []
    for a in patterns: files += glob.glob(os.path.join(temporal_folder,a))
    if len(files)==0:
        qtwrap.show_warning("Nothing to save yet, run a calculation first")
        return
    os.makedirs(save_folder,exist_ok=True) # create the folder if not there
    for f in files: shutil.copy(f,save_folder) # copy this file
    names = ", ".join(sorted(os.path.basename(f) for f in files))
    status("Saved %d files in %s"%(len(files),save_folder))
    qtwrap.show_info("Saved in\n"+save_folder+"\n\n"+names,title="Data saved")


def save_parameters():
    """Write all the fields of the interface to a JSON file chosen by the user"""
    default = os.path.join(original_path,"tranci_parameters.json")
    name,_ = QtWidgets.QFileDialog.getSaveFileName(window,"Save parameters",
              default,"Tranci parameters (*.json)")
    if not name: return # cancelled
    if not name.endswith(".json"): name += ".json"
    write_parameters(name)
    status("Parameters saved in "+name)


def write_parameters(name):
    """Write all the fields of the interface to a JSON file"""
    d = qtwrap.get_all_inputs()
    with open(name,"w") as f: json.dump(d,f,indent=2,sort_keys=True)


def load_parameters():
    """Fill the fields of the interface from a JSON file chosen by the user"""
    name,_ = QtWidgets.QFileDialog.getOpenFileName(window,"Load parameters",
              original_path,"Tranci parameters (*.json);;All files (*)")
    if not name: return # cancelled
    with open(name) as f: d = json.load(f)
    if not isinstance(d,dict):
        raise ValueError(name+" is not a tranci parameters file")
    unknown = qtwrap.set_all_inputs(d)
    status("Parameters loaded from "+name)
    if len(unknown)>0:
        qtwrap.show_warning("These entries were ignored, they do not match "
                "any field of this version:\n"+", ".join(unknown))


restart_tranci() # restart tranci


def open_file(name):
    """Open a file with the default application of the system"""
    if sys.platform=="win32": os.startfile(name) # Windows system
    elif sys.platform=="darwin": subprocess.Popen(["open",name]) # Mac system
    else: subprocess.Popen(["xdg-open",name]) # Linux system


def show_pdf():
    name = os.path.abspath("spectrum_ci.pdf") # full path to the pdf
    if not os.path.isfile(name): # no pdf generated yet
      qtwrap.show_warning("No PDF summary found yet, run the calculation first")
      return
    status("Opening "+name)
    open_file(name)


def show_manual():
    name = os.path.abspath(os.path.join(tranciroot,"doc","tranci_manual.pdf"))
    if not os.path.isfile(name):
      qtwrap.show_warning("The manual was not found at "+name)
      return
    open_file(name)


def run_pdflatex():
  """Compile the latex summary, returning the error text, or "" on success"""
  r = subprocess.run(["pdflatex","-interaction=nonstopmode","spectrum_ci.tex"],
          stdout=subprocess.PIPE,stderr=subprocess.STDOUT) # capture the output
  if r.returncode!=0: # show why, instead of claiming success later
    out = r.stdout.decode("utf-8","replace") if r.stdout else ""
    tail = out.splitlines()[-20:]
    print("pdflatex failed (exit %d). Last lines of its output:"%r.returncode)
    for l in tail: print("  "+l)
    return "\n".join(tail) or ("pdflatex exited with code %d"%r.returncode)
  return ""



# entries of the "Parameter to sweep" combobox: attribute of the parameter
# object and unit of the axis; the angles are in units of pi (see get_b)
SWEEP = {
  "SOC": ("soc","eV"),
  "U": ("U","multiplier"),
  "Uniaxial z^2": ("D","eV"),
  "Shear x^2-y^2": ("E","eV"),
  "Octahedral x^4+y^4+z^4": ("O","eV"),
  "Trigonal ((x+y+z)/sqrt(3))^2": ("trigonal","eV"),
  "Uniaxial' z^4": ("z4","eV"),
  "Shear' (xy)^2+(yx)^2": ("x2y2","eV"),
  "B": ("babs","eV"),
  "Theta_B": ("theta_b","$\\pi$ rad"),
  "Phi_B": ("phi_b","$\\pi$ rad"),
  "J": ("jabs","eV"),
  "Theta_J": ("theta_j","$\\pi$ rad"),
  "Phi_J": ("phi_j","$\\pi$ rad"),
}


def sweep_label(stype):
  """Axis label for a swept variable"""
  return stype+" ["+SWEEP[stype][1]+"]"


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
  update_fields(p)
  return p


def update_fields(p):
  """Recompute the vector fields from modulus and angles"""
  p.b = get_b(p.babs,p.theta_b,p.phi_b) # get the magnetic field
  p.j = get_b(p.jabs,p.theta_j,p.phi_j) # get the exchange field


def initialize_one_shot():
  """ Initialize the one shot calculation"""
  with busy("Running the calculation..."):
    p = read_inputs() # read all the inputs
    write_parameters(parameters_file) # keep the inputs next to the results
    at = get_atom() # read the basis from file
    at.wavefunction_z_axis = [get("zaxis_x"),get("zaxis_y"),get("zaxis_z")]
    m = hamiltonians.build_hamiltonian(at,p) # get the hamiltonian
    if do_check:  check_all(at) # check the hamiltonian
    header = hamiltonians.latex_DE(at,p) # string for the hamiltonian
    ls = lowest_states(m,atom=at) # create the object
    ls.disentangle_manifolds(at.jz) # disentangle manifold
    if qtwrap.is_checked("fit_heff"): nw = int(get("nwf_heff"))
    else: nw = None
    if nw is not None: status("Fitting the effective Hamiltonian...")
    write.write_all(ls,header=header,n=nw) # write Latex file
    if shutil.which("pdflatex") is None: # no latex in this system
      status("Calculation done, spectrum_ci.tex written (no pdflatex found)")
      qtwrap.show_warning("pdflatex was not found, so only spectrum_ci.tex "
        "was written.\nInstall a LaTeX distribution (MiKTeX on Windows, "
        "TeX Live or MacTeX otherwise) to get the PDF summary.")
      return
    status("Compiling the PDF summary...")
    err = run_pdflatex() # compile the latex file
    err2 = run_pdflatex() # do it twice, so that the index is right
    err = err or err2
    if err or not os.path.isfile("spectrum_ci.pdf"): # do not claim success
      status("Calculation done, but the PDF summary could not be created")
      qtwrap.show_warning("spectrum_ci.tex was written but pdflatex failed. "
        "Last lines of its output:\n\n"+err)
      return
  status("Done: PDF summary created, press Show pdf to open it")



def get_atom():
  p = read_inputs() # read all the inputs
  from tranci import atom
  at = atom.get_atom(ne=p.n)
  hamiltonians.tol = np.max([1e-8,get("tol_ene")]) 
  hamiltonians.ntol = -int(round(np.log10(hamiltonians.tol)))
  return at # return atom




def initialize_sweep():
  """Launch a sweeping calculation"""
  p = read_inputs() # read all the inputs
  at = get_atom()
  stype = getbox("sweep_variable") # get the variable
  if stype not in SWEEP: raise ValueError("Unknown sweep variable "+str(stype))
  attr = SWEEP[stype][0]
  def fsweep(x):
    """Function to perform the sweep"""
    setattr(p,attr,x) # set the swept parameter
    update_fields(p) # the vector fields depend on modulus and angles
    m = hamiltonians.build_hamiltonian(at,p) # get the hamiltonian
    ls = lowest_states(m,atom=at) # perform the calculation 
    return ls # return the object
  return fsweep # return function


def run_sweep(fsweep,xs):
  """Evaluate fsweep on every point, reporting progress in the status bar"""
  progress.setRange(0,len(xs))
  progress.setValue(0)
  progress.show()
  out = []
  for i,x in enumerate(xs):
    out.append(fsweep(x))
    progress.setValue(i+1)
    status("Sweep: point %d of %d"%(i+1,len(xs)))
  return out


def sweep_states():
  """Common start of every sweep task: the grid and the states on it"""
  fsweep = initialize_sweep() # get the generator function
  xs = get_sweep_parameters() # get the array
  return xs,run_sweep(fsweep,xs)


def plot_eigenvalues(write=True,center=True):
  """Plots the excited states"""
  ###############################
  ###############################
  with busy("Sweeping..."):
    xs,gst = sweep_states()
  ys = [g.evals_full for g in gst] # get all the eigenvalues
  fig = py.figure() # create figure
  fig.subplots_adjust(.2,.15) # adjust the subplots
  ys = np.array(ys).transpose() # row is same eigenvector evolving
  # number of energies to plot
  nenergies = int(get("nplot")) # number of energies to plot
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
  status("Sweep done, data written to EIGENVALUES.OUT")
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
  with busy("Sweeping..."):
    xs,gst = sweep_states()
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
  status("Sweep done, data written to DEGENERACY.OUT")
  py.show()



def get_operator(at,oname):
  """Many-body operator selected in the Operator combobox"""
  ops = {
    "Sx": at.sx, "Sy": at.sy, "Sz": at.sz,
    "Lx": at.lx, "Ly": at.ly, "Lz": at.lz,
    "Jx": at.jx, "Jy": at.jy, "Jz": at.jz,
    "L2": at.l2, "S2": at.s2, "J2": at.j2, "LS": at.ls,
    "x2": at.x2, "y2": at.y2, "z2": at.z2,
    "x4+y4+z4": at.x4+at.y4+at.z4,
    "S_(111)": (at.sz + at.sx + at.sy)/np.sqrt(3),
    "L_(111)": (at.lz + at.lx + at.ly)/np.sqrt(3),
  }
  if oname not in ops: raise ValueError("Unknown operator "+str(oname))
  return ops[oname]


def plot_operator():
  """Plots the eigenvalues of a certain operator in the GS"""
  ###############################
  ###############################
  with busy("Sweeping..."):
    xs,gst = sweep_states()
    at = get_atom() # get the atom object
    oname = getbox("operator_name")
    op = get_operator(at,oname)
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
  py.ylabel("$\\langle "+oname+"\\rangle$")  # label for the y axis
  stype = getbox("sweep_variable")
  py.xlabel(sweep_label(stype))  # label for the x axis
  fig.set_facecolor("white")
  py.tight_layout()
  status("Sweep done, data written to OPERATOR_VALUES.OUT")
  py.show() # show graph



# Add the tranci logo
from qtwrap import set_logo
tranci_logo = tranciroot+"/logos/orbitals.png"
set_logo("logo",tranci_logo)

# second-quantized formula of every parameter, rendered with mathtext
import formulas
formulas.add_formulas(window)


# create signals: buttons and menu actions
signals = dict()
signals["initialize_one_shot"] = initialize_one_shot  # initialize and run
signals["plot_spectrum"] = plot_spectrum  # initialize and run
signals["plot_excitations"] = plot_excitations  # initialize and run
signals["plot_degeneracy"] = plot_degeneracy  # initialize and run
signals["plot_operator"] = plot_operator  # initialize and run
signals["show_pdf"] = show_pdf  # show pdf with the results
signals["save_tranci"] = save_tranci  # copy the results to ./tranci_data
signals["menu_show_pdf"] = show_pdf
signals["menu_save_tranci"] = save_tranci
signals["load_parameters"] = load_parameters
signals["save_parameters"] = save_parameters
signals["show_manual"] = show_manual
signals["quit_tranci"] = window.close

window.connect_clicks(signals) 
status("Ready. Set the parameters and press Initialize and run.")
window.run()
