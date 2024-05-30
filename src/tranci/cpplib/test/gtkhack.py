
from input import calculation

execfile("interface_ci.py") # execute the interface


if calculation=="one_shot":
  initialize_one_shot(0)
  show_pdf(0)
elif calculation=="plot_excited":
  plot_excited(0)
elif calculation=="plot_degeneracy":
  plot_degeneracy(0)
elif calculation=="plot_operator":
  plot_operator(0)

