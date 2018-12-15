#!/usr/bin/python

# run the calculation for all the occupations

name = "nelectrons.in" # input file

import os

os.system("rm -r cilib") # remove final folder
os.system("mkdir cilib") # create final folder

os.chdir("src/") # go to folder
#for i in range(1,10): # loop over occupation
for i in range(1,10): # loop over occupation
  f = open(name,"w")
  f.write(str(i)+"\n") # write number of electrons
  f.close()
  os.system("./main.x") # run calculation
  os.system("mkdir "+str(i))
  os.system("cp *.op "+str(i))
  os.system("cp *.out "+str(i))
  os.system("cp *.in "+str(i))
  os.system("mv "+str(i)+" ../cilib") # move all the cilib
  
