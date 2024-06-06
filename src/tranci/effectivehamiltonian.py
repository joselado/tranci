import numpy as np
from . import hamiltonians
#from numba import jit

import jax
jax.config.update('jax_platform_name', 'cpu')
from jax import jit
from jax import grad
import jax.numpy as jnp

#@jit(nopython=True)
def errorf(v,diff,ms,simp=1e1,cutoff=1e-6): # function to minimize
    """Error function"""
    n = len(ms) # number of matrices
    rv = v[0:n] # real part
#    iv = v[n:2*n] # imaginary part
#    zv = rv+1j*iv # complex vector
    diff = diff.copy()
    zv = rv
    for i in range(len(ms)): # loop over ms
        diff = diff - zv[i]*ms[i] # add this contribution
    error = jnp.mean(jnp.abs(diff)**2) # error
    zv = zv[1:] # all except identity
    coef = jnp.abs(zv)/jnp.sum(jnp.abs(zv)) # normalize
#    coef = coef[coef>1e-8] # only big enough
#    error = (cutoff+error)*(1.0 - simp*jnp.sum(coef*jnp.log(coef)))
    return error

errorf_jax = jit(errorf)
jacobian_jax = jit(grad(errorf,argnums=0))

def fit_matrix(h,d,cutoff=1e-4,ntries=40,simp = 1e1):
    """Fit a matrix with a dictionary of matrices"""
    ms = np.array([d[key] for key in d]) # redefine as array
    n = len(ms) # number of matrices
    mh = h.copy() # make a copy of the Hamiltonian
    def f(v): # function to minimize
        return errorf(v,mh,ms,simp=simp)
    def jac(v):
        return jacobian_jax(v,mh,ms,simp=simp)
    from scipy.optimize import minimize
    def fopt(): # perform one minimization
        x0 = np.random.random(n)-.5 # random guess
        sol = minimize(f,x0,jac=jac) # with the Jacobian
#        sol = minimize(f,x0,method="Powell",
#                options={'xtol': 1e-6, 'ftol': 1e-6,
#                    'maxiter': 100000,
#                    'maxfev': 100000})
        x = sol.x # solution of the minimization
        x = x[0:n] #+ 1j*x[n:2*n] # redefine as complex
        error = f(x) # compute error
#        print("Error",error)
        return error,x # return solution
    outs = [fopt() for i in range(ntries)] # compute several
    x = [ix for (iy,ix) in sorted(outs,key=lambda x: x[0])][0] # take the smallest one
    error_min = np.min([e for (e,ix) in outs])
    print("Minimum fitting error",error_min)
    errors = []
    out = dict()
    ii = 0
    h0 = np.zeros(mh.shape)
    for key in d: # loop over the operators 
        if np.abs(x[ii])>cutoff:
          out[key] = x[ii]
          h0 = h0 + x[ii]*d[key]
#          print(x[ii])
#          print(np.round(d[key],2))
        ii += 1 # increase counter
    print("Original eigenvalues")
    print(np.round(np.linalg.eigvalsh(h),6))
    print("New eigenvalues")
    print(np.round(np.linalg.eigvalsh(h0),6))
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


def get_s_operators(atom):
    """Return the SJ operators"""
    dd = dict()
    dd["\\hat S_x"] = atom.sx
    dd["\\hat S_y"] = atom.sy
    dd["\\hat S_z"] = atom.sz
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





def get_fitting_operators(lowest,nt=2,n=2,npow=4,dd=None):
    atom = lowest.atom # get the atom object
    if dd is None: dd = get_ls_operators(atom)
    out = dict() # dictionary
    iden = np.identity(atom.lz.shape[0]) # identity
    out[("Id")] = lowest.get_representation(iden,n=n)
    # linear terms
    for ip in range(npow):
      if nt>0: # linear terms
        for di in dd: # loop
          m = dd[di] # store this term
          for ii in range(ip-1): m = m@m # power
          if ip==0: spow = ""
          else: spow = "^"+str(ip+1)
          m = lowest.get_representation(m,n=n)
  #        out[(di)] = m # store this term
          if acceptable_matrix(m,out): # if the matrix can be accepted
            out[(di+spow)] = m.copy() # store this term
    # bilinear terms
      if nt>1: # bilinear terms
        for di in dd: # loop
          for dj in dd: # loop
            mi = dd[di]
            mj = dd[dj]
            for ii in range(ip-1): 
              mi = mi@mi # power
              mj = mj@mj # power
            if ip==0: spow = ""
            else: spow = "^"+str(ip+1)
            m = mi@mj
            m = lowest.get_representation(m,n=n)
  #          out[(di,dj)] = m # store this matrix
            if acceptable_matrix(m,out): # if the matrix can be accepted
              out[(di+spow,dj+spow)] = m.copy() # store this matrix
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
    out += "H = \n"+zform(cmax)+" \\left [ " # output string
    ik = 0 # counter
    nk = len(keys)
    for key in keys: # loop
        c = np.round(d[key]/cmax,4) # round the number
#        print(c,key)
        if np.abs(c)<tol: continue
        if .99<c<1.01: out += "  "
        else: out += zform(c) + "  " # normalize
        out += key2latex(key) # create the name
        ik +=1 # increase counter
        if ik<nk: out += " + \n" # new line
        if ik%3==0: 
          out += "\\\\ \n" # new line
    out += " \\right ] \n" # last line
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



def effective_spin_hamiltonian(lowest,H=None,n=2,nt=2,operators=None):
    """Compute the effective Hamiltonian in Latex form"""
    # get the Hmailtonian
    if H is None: H = lowest.h
    h = lowest.get_representation(H,n=n) # Hamiltonian
#    print(h) ; exit()
    h = h - np.identity(h.shape[0])*np.trace(h)/h.shape[0] # no trace
    atom = lowest.atom # get the atom object
    text = "Hamiltonian written in the low energy manifold with "+str(n)+" states\n"
    if operators is None:
        ops = get_s_operators(atom) # LJ operators
    else:
        ops = dict()
        for key in operators: ops[key] = atom.Operator[key]
    out = get_fitting_operators(lowest,nt=nt,n=n,dd=ops) # get the operators
    # project onto the desired low energy manifold
    # now fit the Hamiltonian
    coef = fit_matrix(h,out) # fit the matrix and return dictionary
    try: del coef[("Id")]
    except: pass
    if len(coef)==0: return ""
    text += "\\begin{equation}\n"
    text +=  dict2latex(coef) # return the latex format
    text += "\\end{equation}\n\n"
    return text
    from .write import matrix2latex
    ops = dict() # dictionary with effective operators
    dd = get_lsj_operators(atom) # get the LSJ operators
    for key in dd: # write all the operators
        m = lowest.get_representation(dd[key],n=n)
        ops[key] = m # save
    from .latexalgebra import effective_algebra
    text += effective_algebra(ops) # write down the effective algebra
    return text


