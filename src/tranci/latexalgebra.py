import numpy as np
from .write import matrix2latex


def effective_algebra(dd):
    """Add the effective algebra of the operators"""
    text = "\\subsection{Projected operators}\n\n"
    text += "Spin, orbital and total angular momentum operators written in "
    text += "the low energy manifold.\n\n"
    for key in dd: # write all the operators
        text += matrix2latex(dd[key],name=key2latex(key)) # get this matrix
    comm = get_commutations(dd)
    if comm=="": return text
    text += "\\subsection{Commutation relations}\n\n"
    text += "Commutators and squares of the projected operators that are "
    text += "proportional to another projected operator.\n\n"
    text += "\\begin{align*}\n"+comm+"\\end{align*}\n\n"
    return text

def return_commutation(keyi,keyj,m,dd):
    out = ""
    for keyk in dd: # loop
        if is_proportional(m,dd[keyk]):
            c = ratio(m,dd[keyk]) # ratio
            out += "  \\big[ "+key2latex(keyi)+", "+key2latex(keyj)+" \\big] &= "
            out += zform(c)+"\\,"+key2latex(keyk)+" \\\\\n"
    return out



def return_square(keyi,m,dd):
    out = ""
    for keyk in dd: # loop
        if is_proportional(m,dd[keyk]):
            c = ratio(m,dd[keyk]) # ratio
            out += "  "+key2latex(keyi)+"^2 &= "
            out += zform(c)+"\\,"+key2latex(keyk)+" \\\\\n"
    return out



def get_commutations(dd):
    """Return the rows of an align* block with every commutator or square
    that is proportional to one of the operators, empty if there is none"""
    out = "" # empty string
    keys = list(dd)
    for (ii,keyi) in enumerate(keys): # loop
      mi = dd[keyi] # this matrix
      for keyj in keys[ii+1:]: # [A,B] = -[B,A], write each pair once
          mj = dd[keyj] # this other matrix
          m = mi@mj - mj@mi # commutator
          out += return_commutation(keyi,keyj,m,dd)
    for keyi in keys: out += return_square(keyi,dd[keyi]@dd[keyi],dd)
    return out


def matrix2vector(m):
    """Convert a matrix into a vector"""
    n = m.shape[0] # get the dimension
    v = np.zeros(n**2,dtype=np.complex128) # to a vector
    v[0:n**2] = m.reshape(n**2)
    return v


def is_proportional(a,b):
    """Check if two vectors are proportional"""
    a = matrix2vector(a)
    b = matrix2vector(b)
    out = braket(a,b)
    aa = np.sqrt(braket(a,a))
    bb = np.sqrt(braket(b,b))
    if np.abs(aa)<1e-6: return False
    if np.abs(bb)<1e-6: return False
    if np.abs(out)<1e-6: return False
    #/(np.sqrt(braket(a,a))*np.sqrt(braket(b,b)))
    out = out/(aa*bb) # normalize
    if 0.9<np.abs(out)<1.1: return True
    return False


def ratio(a,b):
    """Ratio between two vectors"""
    a = matrix2vector(a)
    b = matrix2vector(b)
    out = braket(b,a)/braket(b,b) # coefficient c such that a = c*b
    return out




def braket(a,b):
    """Compute braket"""
    return np.conjugate(a).dot(b)



def key2latex(key):
    if type(key)==str: return key
    return " ".join(key)



from .numberformat import zform
