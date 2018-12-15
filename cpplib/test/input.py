
# Everythong is in eV


##########################
##  Begining of inputs
##########################
D = 0.0 # uniaxial crystal field
E = 0.0 # 
O = 0.0 # octahedral crystal field
U = 7.0 # electron electron interaction
n = 2 # number of electrons
soc = 0.0 # spin orbit coupling
z4 = 0.0 # cuartic axial symmetry
B = 0.0  # type of sweep
theta = 0.0 # theta angle of the magnetic field
phi = 0.0 # phi angle of the magnetic field

# sweep parameters # (see calculation)
# in the sweep all the parameters are kept fixed except the one selected
sweeping_type = "D" # type of sweep
initial = 0.0 # initial value in the sweep
final = 0.4 # final value in the sweep
steps = 20 # number of steps in the sweep
######################


# opertor to plot #
# Valid options  Sx,Sy,Sz,Lx,Ly,Lz,Jx,Jy,Jz,L2,S2,J2
operator = "S2" # name of the operator to plot (see calculation)




#############################################
# choose between the following calculations #
#############################################

# comment all except the one you want

calculation = "one_shot"  # generate pdf for the ground state

# affected by sweep
calculation = "plot_excited" # plot eigenvalues of excited state
calculation = "plot_operator" # plot expectation value over the GS manifold
#calculation = "plot_degeneracy"  # degeneracy of ground state


