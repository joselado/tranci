#!/usr/bin/env python3

# run the calculation for all the occupations

name = "nelectrons.in" # input file

import os
import sys
import glob
import shutil
import subprocess

name_exe = "main.exe" if sys.platform=="win32" else "main.x" # name of the executable
exe = os.path.join(".",name_exe) # executable in the current folder

if os.path.exists("cilib"): shutil.rmtree("cilib") # remove final folder
os.mkdir("cilib") # create final folder

os.chdir("src") # go to folder
if not os.path.isfile(name_exe): # the generator must have been compiled first
  raise SystemExit("%s not found in src/, compile it with ./compile.sh first"%name_exe)
for i in range(1,10): # loop over occupation
  # remove the previous occupation's output, so a partial run cannot be
  # mistaken for a complete one
  for pattern in ["*.op","*.out"]:
    for g in glob.glob(pattern): os.remove(g)
  f = open(name,"w")
  f.write(str(i)+"\n") # write number of electrons
  f.close()
  subprocess.run([exe],check=True) # run calculation, abort if it fails
  folder = os.path.join("..","cilib",str(i)) # final folder for this occupation
  os.mkdir(folder) # create the folder
  for pattern in ["*.op","*.out","*.in"]: # loop over the generated files
    for g in glob.glob(pattern): shutil.copy(g,folder) # copy this file
