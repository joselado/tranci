import numpy as np
from .write import matrix2latex


def effective_algebra(dd):
    """Add the effective algebra of the operators"""
    text = ""
    for key in dd: # write all the operators
        text += matrix2latex(dd[key],name=key) # get this matrix
    text += "\\subsection{Commutations}\n\n" # empty string
    text += get_commutations(dd)
    return text

def return_commutation(keyi,keyj,m,dd):
    out = ""
    for keyk in dd: # loop
        if is_proportional(m,dd[keyk]):
            out += "\\begin{equation}\n"
            out += key2latex(keyi) +"  "
            out += key2latex(keyj) +"  -"
            out += key2latex(keyj) +"  "
            out += key2latex(keyi) +"  ="
            c = ratio(m,dd[keyk]) # ratio
            out += zform(c)
            out += key2latex(keyk) +"  \n"
            out += "\\end{equation}\n"
    return out



def return_square(keyi,m,dd):
    out = ""
    for keyk in dd: # loop
        if is_proportional(m,dd[keyk]):
            out += "\\begin{equation}\n"
            out += key2latex(keyi) +"  "
            out += key2latex(keyi) +"  ="
            c = ratio(m,dd[keyk]) # ratio
            out += zform(c)
            out += key2latex(keyk) +"  \n"
            out += "\\end{equation}\n"
    return out



def get_commutations(dd):
    """Check out if there is any interesting commutation relation"""
    out = "" # empty string
    for keyi in dd: # loop
      mi = dd[keyi] # this matrix
      for keyj in dd: # loop
          mj = dd[keyj] # this other matrix
          m = mi@mj - mj@mi # commutator
          out += return_commutation(keyi,keyj,m,dd)
      out += return_square(keyi,mi@mi,dd)
    return out


def matrix2vector(m):
    """Convert a matrix into a vector"""
    n = m.shape[0] # get the dimension
    v = np.zeros(n**2,dtype=np.complex) # to a vector
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
    print(np.abs(out))
    if 0.9<np.abs(out)<1.1: return True
    return False


def ratio(a,b):
    """Check if two vectors are proportional"""
    a = matrix2vector(a)
    b = matrix2vector(b)
    out = braket(a,b)
    aa = np.sqrt(braket(a,a))
    bb = np.sqrt(braket(b,b))
    out = out/(aa*bb) # normalize
    return out




def braket(a,b):
    """Compute braket"""
    return np.conjugate(a).dot(b)



def key2latex(key):
    if type(key)==str: return key
    out = ""
    for k in key: out += k + "  "
    return out



from .numberformat import zform
