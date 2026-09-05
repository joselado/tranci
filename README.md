# DESCRIPTION #
This program allows to perform CI calculations for d orbitals, by modifying number of electrons crystal fields and spin orbit coupling.

It provides a variaty of results in pdf format, and allows to sweep over the
different parameters of the system.

# DISCLAIMER #
This is still a version under heavy development.

# HOW TO INSTALL #
Execute the script install.py

# USAGE #
Execute "tranci" in a terminal

# DEPENDENCIES #
This library requires Python3, together with the libraries
  - numpy
  - scipy
  - matplotlib
  - PyQt5 (only for the graphical interface)
  - jax (only for the effective Hamiltonian fitting)

The "pdflatex" command should be installed in the system to show the pdf summary.
In Windows this is provided by MiKTeX, in Linux and Mac by TeX Live or MacTeX.

# COMPATIBILITY #
This program works in Linux, Mac and Windows.

In Linux and Mac, execute "python install.py" and afterwards run "tranci" in a
terminal. In Windows, execute "python install.py" and open a new terminal, or
run "bin\tranci.bat" directly.

