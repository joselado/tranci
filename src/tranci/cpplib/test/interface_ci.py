from gi.repository import Gtk as gtk
builder = gtk.Builder()
#from workaround import gtk,builder
import numpy as np
from read import read_matrix
from atom import CIatom
from hamiltonians import DE,eigenvalues,lowest_states,latex_DE
import os # for calling the terminal
import write # for writing in latex 
import matplotlib.pyplot as py
from check import check_all # check that the hamiltonian is right
do_check = True # perform chec of the hamiltonian


def get(name):
  """Get the value of a certain variable"""
  return float(builder.get_object(name).get_text())


def show_pdf(dummy):
  os.system("xdg-open spectrum_ci.pdf")



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
  p.O = get("O")
  p.U = get("U")
  p.soc = get("soc")
  p.theta = get("theta") # theta for magnetic field
  p.phi = get("phi") # phi for Zeeman
  p.babs = get("B") # absolute Zeeman
  p.b = get_b(p.babs,p.theta,p.phi) # get the magnetic field
  return p

def initialize_one_shot(dummy):
  """ Initialize the one shot calculation"""
  p = read_inputs() # read all the inputs
  os.system("cp ../cilib/"+ str(p.n) + "/* ./") # copy input files
#  os.system("cp "+ str(p.n) + "/* ./") # copy input files
  at = CIatom() # create the CI object
  at.read() # read all the matrices
  at.get_basis() # read the basis from file
  m = DE(at,D=p.D,E=p.E,soc=p.soc,U=p.U,O=p.O,z4=p.z4,b=p.b) # get the hamiltonian
  if do_check:  check_all(at) # check the hamiltonian
  header = latex_DE(at,D=p.D,E=p.E,soc=p.soc,U=p.U,O=p.O,z4=p.z4,n=p.n,b=p.b) # string for the hamiltonian
  ls = lowest_states(m) # create the object
  ls.disentangle_manifolds(at.jz)
  write.write_all(at,ls,header=header)
  os.system("pdflatex spectrum_ci.tex")
  os.system("cp spectrum_ci.pdf ../") # copy to the previous folder



def get_atom():
  p = read_inputs() # read all the inputs
#  os.system("cp "+ str(p.n) + "/* ./") # copy input files
  os.system("cp ../cilib/"+ str(p.n) + "/* ./") # copy input files
  at = CIatom() # create the CI object
  at.read() # read all the matrices
  at.get_basis() # read the basis from file
  return at # return atom




def initialize_sweep():
  """Launch a sweeping calculation"""
  p = read_inputs() # read all the inputs
#  os.system("cp "+ str(p.n) + "/* ./") # copy input files
  os.system("cp ../cilib/"+ str(p.n) + "/* ./") # copy input files
  at = get_atom() # get the atom
  def fsweep(x):
    """Function to perform the sweep"""
    stype = builder.get_object("sweeping_type").get_active_text()
    if stype == "U": m = DE(at,D=p.D,E=p.E,soc=p.soc,U=x,O=p.O,z4=p.z4,b=p.b)
    elif stype == "D": m = DE(at,D=x,E=p.E,soc=p.soc,U=p.U,O=p.O,z4=p.z4,b=p.b)
    elif stype == "E": m = DE(at,D=p.D,E=x,soc=p.soc,U=p.U,O=p.O,z4=p.z4,b=p.b)
    elif stype == "soc": m = DE(at,D=p.D,E=p.E,soc=x,U=p.U,O=p.O,z4=p.z4,b=p.b)
    elif stype == "U": m = DE(at,D=p.D,E=p.E,soc=p.soc,U=x,O=p.O,z4=p.z4,b=p.b)
    elif stype == "z4": m = DE(at,D=p.D,E=p.E,soc=p.soc,U=p.U,O=p.O,z4=x,b=p.b)
    elif stype == "O": m = DE(at,D=p.D,E=p.E,soc=p.soc,U=p.U,O=x,z4=p.z4,b=p.b)
    elif stype == "B": m = DE(at,D=p.D,E=p.E,soc=p.soc,U=p.U,O=p.O,z4=p.z4,b=get_b(x,p.theta,p.phi))
    elif stype == "Theta": m = DE(at,D=p.D,E=p.E,soc=p.soc,U=p.U,O=p.O,z4=p.z4,b=get_b(p.babs,x,p.phi))
    elif stype == "Phi": m = DE(at,D=p.D,E=p.E,soc=p.soc,U=p.U,O=p.O,z4=p.z4,b=get_b(p.babs,p.theta,x))
    else: raise # raise error
    ls = lowest_states(m) # perform the calculation 
    return ls # return the object
  return fsweep # return function

def plot_excited(dummy):
  """Plots the excited states"""
  ###############################
  ###############################
  fsweep = initialize_sweep() # get the generator function
  xs = np.linspace(get("initial"),get("final"),get("steps")) # array for sweep
  gst = [fsweep(ix) for ix in xs]  # create the list of objects
  ys = [g.get_excitations() for g in gst] # get the energies of the excitations
  fig = py.figure() # create figure
#  fig.subplot_adjust()
  for (x,y) in zip(xs,ys):
    py.scatter([x for iy in y],y,c="black") 
  py.xlim([min(xs),max(xs)]) # 
  py.ylabel("Energy [eV]")  # label for the y axis
  stype = builder.get_object("sweeping_type").get_active_text()
  py.xlabel(stype+" [eV]")  # label for the x axis
  fig.set_facecolor("white")
  py.show()





def plot_degeneracy(dummy):
  """Plots the degeneracy of the ground state"""
  ###############################
  ###############################
  fsweep = initialize_sweep() # get the generator function
  xs = np.linspace(get("initial"),get("final"),get("steps")) # array for sweep
  gst = [fsweep(ix) for ix in xs]  # create the list of objects
  ds = [g.get_gs_multiplicity() for g in gst] # get degeneracies 
  fig = py.figure() # create figure
  py.plot(xs,ds,c="green",linewidth=5.0) 
  py.xlim([min(xs),max(xs)]) # 
  py.ylim([0,max(ds)+1]) # 
  py.ylabel("Degeneracy")  # label for the y axis
  stype = builder.get_object("sweeping_type").get_active_text()
  py.xlabel(stype + "[eV]")  # label for the y axis
  fig.set_facecolor("white")
  py.show()





def plot_operator(dummy):
  """Plots the eigenvalues of a certain operator in the GS"""
  ###############################
  ###############################
  fsweep = initialize_sweep() # get the generator function
  xs = np.linspace(get("initial"),get("final"),get("steps")) # array for sweep
  gst = [fsweep(ix) for ix in xs]  # create the list of objects
  at = get_atom() # get the atom object

  ########################
  ########################
  ########################
  oname = builder.get_object("operator").get_active_text()
  if (oname=="Sx"): op = at.sx     # get this operator 
  elif (oname=="Sy"): op = at.sy   # get this operator 
  elif (oname=="Sz"): op = at.sz   # get this operator 
  elif (oname=="Jx"): op = at.jx   # get this operator
  elif (oname=="Jy"): op = at.jy   # get this operator
  elif (oname=="Jz"): op = at.jz   # get this operator
  elif (oname=="Lx"): op = at.lx   # get this operator 
  elif (oname=="Ly"): op = at.ly   # get this operator 
  elif (oname=="Lz"): op = at.lz   # get this operator 
  elif (oname=="L2"): op = at.l2   # get this operator 
  elif (oname=="S2"): op = at.s2   # get this operator 
  elif (oname=="J2"): op = at.j2   # get this operator 
  else: print oname ; raise
  ########################
  ########################
  ########################

  evals = [g.get_gs_projected_eigenvalues(op) for g in gst] # get op eigen 
  fig = py.figure() # create figure
  for (x,y) in zip(xs,evals):
    py.scatter([x for iy in y],y,c="black") 
  py.xlim([min(xs),max(xs)]) # 
#  py.ylim([min(evals),max(evals)]) # 
  py.ylabel(oname)  # label for the y axis
  stype = builder.get_object("sweeping_type").get_active_text()
  py.xlabel(stype+"  [eV]")  # label for the y axis
  fig.set_facecolor("white")
  py.show() # show graph

















# create signals
signals = dict()
signals["on_window_destroy"] = gtk.main_quit  # close the window
signals["initialize_one_shot"] = initialize_one_shot  # initialize and run
signals["initialize_sweep"] = plot_excited  # initialize and run
signals["plot_degeneracy"] = plot_degeneracy  # initialize and run
signals["plot_operator"] = plot_operator  # initialize and run
signals["show_pdf"] = show_pdf  # show pdf with the results

class CIApp(object):
        def __init__(self):
            builder.add_from_file("interface.xml")
            builder.connect_signals(signals)
            self.window = builder.get_object("ci_window")
            self.window.show()



if __name__ == "__main__":
        app = CIApp()
        gtk.main()


