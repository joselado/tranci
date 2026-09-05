#!/usr/bin/env python3

import os
import platform

bindir = os.path.join(os.path.dirname(os.path.realpath(__file__)),"bin") # path to bin


def install_unix():
  """Add tranci to the PATH in Linux and Mac"""
  addpath = "\n\n# Path to tranci, CI for transition metals on surfaces\n"
  addpath += "export PATH=$PATH:"+bindir+"\n\n" # folder with the executable
  if platform.system()=="Linux":
    rcfile = os.path.join(os.path.expanduser("~"),".bashrc") # path to .bashrc
    print("Detected Linux system")
  else:
    rcfile = os.path.join(os.path.expanduser("~"),".bash_profile") # path to .bash_profile
    print("Detected Mac system")
  with open(rcfile,"a") as f: f.write(addpath) # add the line
  print("Added to your $PATH in",rcfile)


def install_windows():
  """Add tranci to the PATH of the user in Windows"""
  import winreg
  print("Detected Windows system")
  key = winreg.OpenKey(winreg.HKEY_CURRENT_USER,"Environment",0,
          winreg.KEY_READ | winreg.KEY_WRITE) # open the user environment
  try: (old,vtype) = winreg.QueryValueEx(key,"Path") # current PATH of the user
  except FileNotFoundError: (old,vtype) = ("",winreg.REG_EXPAND_SZ) # no PATH yet
  folders = [o for o in old.split(";") if o!=""] # folders in the PATH
  if bindir.lower() in [o.lower() for o in folders]: # already there
    print("tranci is already in your PATH")
  else:
    folders.append(bindir) # add the folder with the executable
    winreg.SetValueEx(key,"Path",0,vtype,";".join(folders)) # write the new PATH
    print("Added to your PATH")
  winreg.CloseKey(key) # close the key
  broadcast_path_change() # tell the already open windows about the new PATH
  print("Open a new terminal to use it")


def broadcast_path_change():
  """Tell Windows that the environment changed, otherwise new terminals
  started from Explorer keep using the old cached PATH"""
  import ctypes
  HWND_BROADCAST = 0xFFFF # all the top level windows
  WM_SETTINGCHANGE = 0x1A # message for a change in the settings
  SMTO_ABORTIFHUNG = 0x0002 # do not block if a window is hung
  ctypes.windll.user32.SendMessageTimeoutW(HWND_BROADCAST,WM_SETTINGCHANGE,0,
          "Environment",SMTO_ABORTIFHUNG,5000,None) # send the message


if platform.system()=="Windows": install_windows()
else: install_unix()

print("Now you can execute anywhere tranci")
