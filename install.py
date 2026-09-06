#!/usr/bin/env python3

import os
import platform

bindir = os.path.join(os.path.dirname(os.path.realpath(__file__)),"bin") # path to bin


def get_rcfile():
  """Return the startup file of the login shell, and its name"""
  shell = os.path.basename(os.environ.get("SHELL","")) # login shell
  home = os.path.expanduser("~")
  if shell=="zsh": return os.path.join(home,".zshrc"),"zsh"
  if shell in ("bash",""): # bash, or unknown
    if platform.system()=="Linux": return os.path.join(home,".bashrc"),"bash"
    return os.path.join(home,".bash_profile"),"bash" # login bash on Mac
  # some other shell: fall back to the platform default for bash
  if platform.system()=="Linux": return os.path.join(home,".bashrc"),shell
  return os.path.join(home,".bash_profile"),shell


def install_unix():
  """Add tranci to the PATH in Linux and Mac"""
  # quote the path, otherwise a folder with a space breaks every new shell
  exportline = 'export PATH="$PATH:%s"' % bindir
  addpath = "\n\n# Path to tranci, CI for transition metals on surfaces\n"
  addpath += exportline + "\n\n" # folder with the executable
  print("Detected",platform.system(),"system")
  rcfile,shell = get_rcfile() # startup file of the login shell
  if os.path.isfile(rcfile): # do not append the same line twice
    with open(rcfile,"r") as f: current = f.read()
    if exportline in current:
      print("tranci is already in your $PATH in",rcfile)
      return
  with open(rcfile,"a") as f: f.write(addpath) # add the line
  print("Added to your $PATH in",rcfile,"(startup file of %s)"%shell)


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
