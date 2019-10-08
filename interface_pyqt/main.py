import os
import sys
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

## random name for the temporal folder
temporal_folder = "/tmp/tranci_tmp_"+str(np.random.randint(10000))

# go to a temporal folder
def restart_tranci():
  """Fully restart tranci"""
  tmpfol = temporal_folder
  print("Temporal Tranci folder is",tmpfol)
  os.system("rm -rf "+tmpfol) # remove temporal folder
  os.system("mkdir "+tmpfol) # remove temporal folder
  os.chdir(tmpfol) # go to the temporal folder

def save_tranci():
    """Save the generated data in a file"""
    save_folder = original_path+"/tranci_data" # name of the folder
    def tcopy(a): 
        os.system("cp "+temporal_folder+"/"+a+" "+save_folder+" 2>/dev/null")
    os.system("mkdir "+save_folder+" 2>/dev/null")
    tcopy("*.tex")
    tcopy("*.pdf")
    tcopy("*.OUT")
    print("Saved tranci data in ",save_folder)



restart_tranci() # restart tranci


def show_pdf(dummy):
    import platform
    if platform.system()=="Linux":
      os.system("xdg-open spectrum_ci.pdf")
    else: # Mac system
      os.system("open spectrum_ci.pdf")



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

def initialize_one_shot(dummy):
  """ Initialize the one shot calculation"""
  p = read_inputs() # read all the inputs
#  os.system("cp "+ str(p.n) + "/* ./") # copy input files
  at = get_atom() # read the basis from file
  m = hamiltonians.build_hamiltonian(at,p) # get the hamiltonian
  if do_check:  check_all(at) # check the hamiltonian
  header = hamiltonians.latex_DE(at,p) # string for the hamiltonian
  ls = lowest_states(m,atom=at) # create the object
  ls.disentangle_manifolds(at.jz) # disentangle manifold
#  ls.get_gtensor() # compute gfactor
#  exit()
  write.write_all(ls,header=header) # write Latex file
  os.system("pdflatex spectrum_ci.tex")
  os.system("pdflatex spectrum_ci.tex") # do it twice
  os.system("cp spectrum_ci.pdf ../") # copy to the previous folder



def get_atom():
  p = read_inputs() # read all the inputs
#  os.system("cp "+ str(p.n) + "/* ./") # copy input files
  os.system("cp "+tranciroot+"cilib/"+ str(p.n) + "/* ./") # copy input files
  at = CIatom() # create the CI object
  at.read() # read all the matrices
  at.get_basis() # read the basis from file
  hamiltonians.tol = np.max([1e-8,get("tol_ene")]) 
  hamiltonians.ntol = -int(round(np.log10(hamiltonians.tol)))
  return at # return atom




def initialize_sweep():
  """Launch a sweeping calculation"""
  p = read_inputs() # read all the inputs
#  os.system("cp "+ str(p.n) + "/* ./") # copy input files
  os.system("cp "+tranciroot+"cilib/"+ str(p.n) + "/* ./") # copy input files
  at = get_atom() # get the atom
  def fsweep(x):
    """Function to perform the sweep"""
    stype = getbox("sweep_variable") # get the variable
    if stype == "U": p.U = x
    elif stype == "z^2": p.D = x
    elif stype == "x^2-y^2": p.E = x
    elif stype == "soc": p.soc = x
    elif stype == "z^4": p.z4 = x
    elif stype == "x^2y^2": p.x4y2 = x
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
  nenergies = int(get("num_ene_plot")) # number of energies to plot
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
  if stype in ["Phi","Theta"]: py.xlabel(stype+" [rad]")  # label for the x axis
  else: py.xlabel(stype+" [eV]")  # label for the x axis
  fig.set_facecolor("white")
  py.show()




def plot_spectrum(dummy):
  """ Plot the different eigenvlaues, centered"""
  plot_eigenvalues(center=True)


def plot_excitations(dummy):
  """ Plot the different eigenvlues, shifted to the GS"""
  plot_eigenvalues(center=False)

def get_sweep_parameters():
  return np.linspace(get("initial_value"),get("final_value"),int(get("steps"))) 


def plot_degeneracy(dummy):
  """Plots the degeneracy of the ground state"""
  ###############################
  ###############################
  xs = get_sweep_parameters() # get the array
  fsweep = initialize_sweep() # get the generator function
  gst = [fsweep(ix) for ix in xs]  # create the list of objects
  ds = [g.get_gs_multiplicity() for g in gst] # get degeneracies 
  fig = py.figure() # create figure
  fig.subplots_adjust(.2,.15) # adjust the subplots
  py.plot(xs,ds,c="green",marker="o") 
  np.savetxt("DEGENERACY.OUT",np.array([xs,ds]).T) # save data
  py.xlim([min(xs),max(xs)]) # 
  py.ylim([0,max(ds)+1]) # 
  py.ylabel("Degeneracy")  # label for the y axis
  stype = getbox("sweep_variable")
  if stype in ["Phi","Theta"]: py.xlabel(stype+" [rad]")  # label for the x axis
  else: py.xlabel(stype + "[eV]")  # label for the y axis
  fig.set_facecolor("white")
  py.show()





def plot_operator(dummy):
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
  if stype in ["Phi","Theta"]: py.xlabel(stype+" [rad]")  # label for the x axis
  else: py.xlabel(stype+"  [eV]")  # label for the y axis
  fig.set_facecolor("white")
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

