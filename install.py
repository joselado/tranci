#!/usr/bin/python

import os

addpath = "\n\n# Path to tranci, CI for transition metals on surfaces\n"
addpath += "export PATH=$PATH:"
addpath += os.getcwd()+"/bin\n\n" # current direction 

import platform
if platform.system()=="Linux":
  bashrc = os.environ["HOME"]+"/.bashrc" # path to .bashrc
  print("Detected Linux system")
else:
  bashrc = os.environ["HOME"]+"/.bash_profile" # path to .bashrc
  print("Detected Mac system")

os.system("echo '"+addpath+"' >> "+bashrc)

print("Added to your $PATH, now you can execute anywhere tranci")


