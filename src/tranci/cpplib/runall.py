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
for i in range(1,10): # loop over occupation
  f = open(name,"w")
  f.write(str(i)+"\n") # write number of electrons
  f.close()
  subprocess.run([exe]) # run calculation
  folder = os.path.join("..","cilib",str(i)) # final folder for this occupation
  os.mkdir(folder) # create the folder
  for pattern in ["*.op","*.out","*.in"]: # loop over the generated files
    for g in glob.glob(pattern): shutil.copy(g,folder) # copy this file
