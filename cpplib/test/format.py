# library with different format stuff

tol = 0.001
guess_number = True

import numpy as np

def fform(x):
  """ Floating format"""
  return "{:10.2f}".format(x)


def zform(z):
  """ Format a complex number """
  if not guess_number:
    return fform(z.real)+"+"+fform(z.imag)+"i"
 
  strn = "" # initialice string for the number
  if np.abs(z)<tol: strn = "0" # zero number
  elif np.abs(z.real)<tol: # zero number
    if z.imag<0: strn +="-"+ recognise_number(np.abs(z.imag))+"i" # positive
    else: strn += recognise_number(np.abs(z.imag))+"i" # positive
  else:  # real part
    strn += recognise_number(z.real)
    if np.abs(z.imag)>tol: # if has imaginary part
      if z.imag<0: strn += recognise_number(np.abs(z.imag))+"i" # positive
      else: strn += "-"+recognise_number(np.abs(z.imag))+"i" # positive
  return strn
#  elif z.imag>0: return recognise_number(z.real)+"+i"+recognise_number(z.imag) # positive
#  elif z.imag<0: return recognise_number(z.real)+"-i"+recognise_number(np.abs(z.imag)) # positive


def recognise_number(x):
  """Tries to recognise a real number"""
  if not guess_number:  return "{:10.2f}".format(x)
  xa = np.abs(x) # absolute value
  s2 = 1./np.sqrt(2.)
  sq2 = np.sqrt(2.)
  sq6 = np.sqrt(6.)
  sq23 = np.sqrt(3./2.)
  s3 = 1./np.sqrt(3.)
  if x>0.: strx = ""
  else: strx = "-"
  if np.abs(xa-s2)<tol: strx += "\\frac{1}{\\sqrt{2}}"
  elif np.abs(xa-s3)<tol: strx += "\\frac{1}{\\sqrt{3}}"
  elif np.abs(xa-s2/2.)<tol: strx += "\\frac{1}{2\\sqrt{2}}"
  elif np.abs(xa-1./2.)<tol: strx += "\\frac{1}{2}"
  elif np.abs(xa-3./2.)<tol: strx += "\\frac{3}{2}"
  elif np.abs(xa-5./2.)<tol: strx += "\\frac{5}{2}"
  elif np.abs(xa-5./4.)<tol: strx += "\\frac{5}{4}"
  elif np.abs(xa-2./5.)<tol: strx += "\\frac{2}{5}"
  elif np.abs(xa-6./5.)<tol: strx += "\\frac{6}{5}"
  elif np.abs(xa-7./2.)<tol: strx += "\\frac{7}{2}"
  elif np.abs(xa-1./3.)<tol: strx += "\\frac{1}{3}"
  elif np.abs(xa-1./4.)<tol: strx += "\\frac{1}{4}"
  elif np.abs(xa-1./5.)<tol: strx += "\\frac{1}{5}"
  elif np.abs(xa-sq2)<tol: strx += "\\sqrt{2}"
  elif np.abs(xa-sq6)<tol: strx += "\\sqrt{6}"
  elif np.abs(xa-sq23)<tol: strx += "\\sqrt{\\frac{3}{2}}"
  elif np.abs(xa-int(round(xa,0)))<tol: strx += str(int(round(xa,0)))
  else: strx += "{:10.2f}".format(xa)
  return strx


