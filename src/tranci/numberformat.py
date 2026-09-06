# library with different format stuff

tol = 0.001 # below this a number is printed as zero
guess_tol = 1e-4 # tolerance for recognising closed forms such as 1/sqrt(2)
guess_number = True

import numpy as np

def fform(x,n=3):
  """ Floating format (fixed width, used for the .OUT data files)"""
  return ("{:10."+str(n)+"f}").format(x)


def texnum(x,n=3):
  """Fixed-point number for a LaTeX table: no padding, and a value that
  rounds to zero is printed as 0.000 instead of -0.000"""
  s = ("{:."+str(n)+"f}").format(float(np.real(x)))
  if float(s)==0.0: s = s.replace("-","")
  return s


def zform(z,guess_tol=None):
  """ Format a complex number. guess_tol overrides the tolerance used to
  recognise closed forms such as 1/sqrt(2) """
  if not guess_number:
    return fform(z.real)+"+"+fform(z.imag)+"i"
 
  strn = "" # initialice string for the number
  if np.abs(z)<tol: strn = "0" # zero number
  elif np.abs(z.real)<tol: # zero number
    if z.imag<0: strn +="-"+ recognise_number(np.abs(z.imag),guess_tol)+"i" # positive
    else: strn += recognise_number(np.abs(z.imag),guess_tol)+"i" # positive
  else:  # real part
    strn += recognise_number(z.real,guess_tol)
    if np.abs(z.imag)>tol: # if has imaginary part
      if z.imag<0: strn += "-"+recognise_number(np.abs(z.imag),guess_tol)+"i" # negative
      else: strn += "+"+recognise_number(np.abs(z.imag),guess_tol)+"i" # positive
  return strn
#  elif z.imag>0: return recognise_number(z.real,guess_tol)+"+i"+recognise_number(z.imag) # positive
#  elif z.imag<0: return recognise_number(z.real,guess_tol)+"-i"+recognise_number(np.abs(z.imag),guess_tol) # positive


def recognise_number(x,guess_tol=None):
  """Tries to recognise a real number as a closed form, otherwise prints it
  with two decimals. guess_tol defaults to the module-level guess_tol, and
  never exceeds tol so that a tighter global tolerance still applies"""
  if not guess_number:  return "{:.2f}".format(x)
  xa = np.abs(x) # absolute value
  s2 = 1./np.sqrt(2.)
  sq2 = np.sqrt(2.)
  sq6 = np.sqrt(6.)
  sq23 = np.sqrt(3./2.)
  s3 = 1./np.sqrt(3.)
  if x>=0.: strx = ""
  else: strx = "-"
  if guess_tol is None: guess_tol = globals()["guess_tol"]
  tol = min(guess_tol,globals()["tol"]) # closed forms only when essentially exact
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
  else: strx += "{:.2f}".format(xa)
  return strx


