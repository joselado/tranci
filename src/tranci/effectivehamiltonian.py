import numpy as np
from . import hamiltonians
from numba import jit

@jit(nopython=True)
def errorf(v,diff,ms): # function to minimize
    """Error fucntion"""
    n = len(ms) # number of matrices
    rv = v[0:n] # real part
#    iv = v[n:2*n] # imaginary part
#    zv = rv+1j*iv # complex vector
    zv = rv
    for i in range(len(ms)): # loop over ms
        diff = diff - zv[i]*ms[i] # add this contribution
    error = np.mean(np.abs(diff)**2) # error
    return error



def fit_matrix(h,d,cutoff=1e-4,ntries=10):
    """Fit a matrix with a dictionary of matrices"""
    ms = np.array([d[key] for key in d]) # redefine as array
    n = len(ms) # number of matrices
    def f(v): # function to minimize
        return errorf(v,h.copy(),ms)
    from scipy.optimize import minimize
    def fopt(): # perform one minimization
        x0 = np.random.random(n)-.5 # random guess
        sol = minimize(f,x0,method="Powell",
                options={'xtol': 1e-6, 'ftol': 1e-6,
                    'maxiter': 100000,
                    'maxfev': 100000})
        x = sol.x # solution of the minimization
        x = x[0:n] #+ 1j*x[n:2*n] # redefine as complex
        error = f(x) # compute error
        print("Error",error)
        return error,x # return solution
    outs = [fopt() for i in range(ntries)] # compute several
    x = [ix for (iy,ix) in sorted(outs)][0] # take the smallest one
    errors = []
    out = dict()
    ii = 0
    h0 = 0.0
    for key in d: # loop over the operators 
        if np.abs(x[ii])>cutoff:
          out[key] = x[ii]
          h0 = h0 + x[ii]*d[key]
        ii += 1 # increase counter
    print(np.linalg.eigvalsh(h),"Original Hamiltonian")
#    h0 = h0 + np.conjugate(h0.T)
#    print(np.linalg.eigvalsh(h0/2.),"Computed Hamiltonian")
    return out # return the coefficients

def get_ls_operators(atom):
    """Return the LS operators"""
    dd = dict()
    dd["\\bar S_x"] = atom.sx
    dd["\\bar S_y"] = atom.sy
    dd["\\bar S_z"] = atom.sz
    dd["\\bar L_x"] = atom.lx
    dd["\\bar L_y"] = atom.ly
    dd["\\bar L_z"] = atom.lz
    for d in dd: dd[d] = dd[d].todense()
    return dd


def get_sj_operators(atom):
    """Return the SJ operators"""
    dd = dict()
    dd["\\bar S_x"] = atom.sx
    dd["\\bar S_y"] = atom.sy
    dd["\\bar S_z"] = atom.sz
    dd["\\bar J_x"] = atom.jx
    dd["\\bar J_y"] = atom.jy
    dd["\\bar J_z"] = atom.jz
    for d in dd: dd[d] = dd[d].todense()
    return dd



def get_lj_operators(atom):
    """Return the SJ operators"""
    dd = dict()
    dd["\\bar L_x"] = atom.lx
    dd["\\bar L_y"] = atom.ly
    dd["\\bar L_z"] = atom.lz
    dd["\\bar J_x"] = atom.jx
    dd["\\bar J_y"] = atom.jy
    dd["\\bar J_z"] = atom.jz
    for d in dd: dd[d] = dd[d].todense()
    return dd



def get_lsj_operators(atom):
    """Return the SJ operators"""
    dd = dict()
    dd["\\bar L_x"] = atom.lx
    dd["\\bar L_y"] = atom.ly
    dd["\\bar L_z"] = atom.lz
    dd["\\bar S_x"] = atom.sx
    dd["\\bar S_y"] = atom.sy
    dd["\\bar S_z"] = atom.sz
    dd["\\bar J_x"] = atom.jx
    dd["\\bar J_y"] = atom.jy
    dd["\\bar J_z"] = atom.jz
    for d in dd: dd[d] = dd[d].todense()
    return dd





def get_fitting_operators(lowest,nt=2,n=2,dd=None):
    atom = lowest.atom # get the atom object
    if dd is None: dd = get_ls_operators(atom)
    out = dict() # dictionary
    iden = np.identity(atom.lz.shape[0]) # identity
    out[("Id")] = lowest.get_representation(iden,n=n)
    # linear terms
    if nt>0: # bilinear terms
      for di in dd: # loop
        m = dd[di] # store this term
        m = lowest.get_representation(m,n=n)
        if acceptable_matrix(m,out): # if the matrix can be accepted
          out[(di)] = m # store this term
    # bilinear terms
    if nt>1: # bilinear terms
      for di in dd: # loop
        for dj in dd: # loop
            m = dd[di]@dd[dj]
            m = lowest.get_representation(m,n=n)
            if acceptable_matrix(m,out): # if the matrix can be accepted
              out[(di,dj)] = m # store this matrix
    return out




def effective_hamiltonian(lowest,n=2,nt=2):
    """Compute the effective Hamiltonian in Latex form"""
    # get the Hmailtonian
    h = lowest.get_representation(lowest.h,n=n) # Hamiltonian
    h = h - np.identity(h.shape[0])*np.trace(h)/h.shape[0] # no trace
    atom = lowest.atom # get the atom object
    ls = get_ls_operators(atom) # LS operators
    sj = get_sj_operators(atom) # SJ operators
    lj = get_lj_operators(atom) # LJ operators
    text = "\\section{Effective Hamiltonian}\n\n\n" #
    text += "This is the Hamiltonian written in the low energy manifold with "+str(n)+" states\n"
    ops = [ls,sj,lj] # operators
    names = ["LS","SJ","LJ"] # names
    for (dd,name) in zip(ops,names): # loop over pairs of effective operators
      out = get_fitting_operators(lowest,nt=nt,n=n,dd=dd) # get the operators
      # project onto the desired low energy manifold
      # now fit the Hamiltonian
      coef = fit_matrix(h,out) # fit the matrix and return dictionary
      try: del coef[("Id")]
      except: pass
      if len(coef)==0: return ""
      text += "\\subsection{Low energy Hamiltonian with "+name+" operators}"
      text += "\\begin{equation}\n"
      text +=  dict2latex(coef) # return the latex format
      text += "\\end{equation}\n\n"
      from .write import matrix2latex
      ops = dict() # dictionary with effective operators
    dd = get_lsj_operators(atom) # get the LSJ operators
    for key in dd: # write all the operators
        m = lowest.get_representation(dd[key],n=n)
        ops[key] = m # save
    from .latexalgebra import effective_algebra
    text += effective_algebra(ops) # write down the effective algebra
    return text



def key2latex(key):
    out = ""
    for k in key: out += k + "  "
    return out


def dict2latex(d,tol=1e-2):
    """Transform the dictionary into a latex form"""
    cs = [d[key] for key in d] # coefficients
    cmax = [iy for (ix,iy) in sorted(zip(np.abs(cs),cs))][-1] 
    keys = [key for key in d] # get the keys
    keys = [iy for (ix,iy) in sorted(zip(-np.abs(cs),keys))] # sort the keys
    out = "\\begin{aligned}\n"
    out += "H = \n"+zform(cmax)+" [ \\\\ \n" # output string
    ik = 0 # counter
    for key in keys: # loop
        c = np.round(d[key]/cmax,3) # round the number
        print(c,key)
        if np.abs(c)<tol: continue
        if .99<c<1.01: out += "  "
        else: out += zform(c) + "  " # normalize
        out += key2latex(key) # create the name
        out += " + \n" # new line
        ik +=1 # increase counter
        if ik%3==0: 
          out += "\\\\ \n" # new line
    out += " ] \n" # last line
    out += "\\end{aligned}\n"
    return out


from .numberformat import zform


def acceptable_matrix(m,ops):
    """Check if it is ok to keep this matrix"""
    if np.sum(np.abs(m))<1e-7: return False
    v = matrix2vector(m)
    out = [v] # list
    for key in ops: # loop over the other matrices
        o = ops[key] # get the matrix
        vo = matrix2vector(o) # convert to vector
        out.append(vo)
    r = np.linalg.matrix_rank(np.array(out),tol=1e-3)
    if r==(len(ops)+1): return True
    else: False
#        proj = braket(v,vo)/(np.sqrt(braket(v,v))*np.sqrt(braket(vo,vo)))
#        if np.abs(proj)>0.98:
      #      print("Skipping")
#            return False
#    return True

def braket(a,b):
    return np.abs(np.conjugate(a).dot(b))

from .latexalgebra import matrix2vector
