import numpy as np
from . import hamiltonians
#from numba import jit

try: # jax is only needed to fit the effective Hamiltonian
    import jax
    jax.config.update('jax_platform_name', 'cpu')
    from jax import jit
    from jax import grad
    import jax.numpy as jnp
except ImportError:
    raise ImportError("Fitting the effective Hamiltonian requires jax, "
            "install it with 'pip install jax'")

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

def fit_matrix(h,d,cutoff=1e-5,ntries=40,simp = 1e1):
    """Fit a matrix with a dictionary of matrices"""
    ms = np.array([d[key] for key in d]) # redefine as array
    n = len(ms) # number of matrices
    mh = h.copy() # make a copy of the Hamiltonian
    def f(v): # function to minimize
        return errorf_jax(v,mh,ms,simp=simp)
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
    e_in = np.linalg.eigvalsh(h) # spectrum to reproduce
    e_out = np.linalg.eigvalsh(h0) # spectrum of the fitted Hamiltonian
    print("Original eigenvalues")
    print(np.round(e_in,6))
    print("New eigenvalues")
    print(np.round(e_out,6))
    scale = np.max(np.abs(e_in)) # energy scale of the manifold
    if scale>0.: # relative error of the fitted spectrum
        fit_matrix.relative_error = float(np.max(np.abs(e_in-e_out))/scale)
    else: fit_matrix.relative_error = 0.
    if fit_matrix.relative_error>1e-2: # the basis cannot represent this H
        print("WARNING: the effective Hamiltonian reproduces the spectrum only "
              "to %.1f%%"%(100*fit_matrix.relative_error))
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





def get_fitting_operators(lowest,nt=2,n=2,npow=2,dd=None):
    atom = lowest.atom # get the atom object
    if dd is None: dd = get_ls_operators(atom)
    out = dict() # dictionary
    iden = np.identity(atom.lz.shape[0]) # identity
    out[("Id")] = lowest.get_representation(iden,n=n)
    # linear terms
    for ip in range(npow):
      if nt>0: # linear terms
        for di in dd: # loop
          m = np.linalg.matrix_power(dd[di],ip+1) # power ip+1
          if ip==0: spow = ""
          else: spow = "^"+str(ip+1)
  #        m = lowest.get_representation(m,n=n)
  #        out[(di)] = m # store this term
          if acceptable_matrix(m,out): # if the matrix can be accepted
            out[(di+spow)] = m.copy() # store this term
    # bilinear terms
      if nt>1: # bilinear terms
        for di in dd: # loop
          for dj in dd: # loop
            mi = np.linalg.matrix_power(dd[di],ip+1) # power ip+1
            mj = np.linalg.matrix_power(dd[dj],ip+1) # power ip+1
            if ip==0: spow = ""
            else: spow = "^"+str(ip+1)
            # symmetrize: for non-commuting operators mi@mj is not Hermitian,
            # so a real-coefficient fit could not represent it and the emitted
            # formula was not Hermitian either
            m = (mi@mj + mj@mi)/2.
  #          m = lowest.get_representation(m,n=n)
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
      # project onto the low energy manifold first: get_fitting_operators mixes
      # these with an n x n identity, so full-size operators cannot be used
      ddp = dict()
      for key in dd: ddp[key] = lowest.get_representation(dd[key],n=n)
      out = get_fitting_operators(lowest,nt=nt,n=n,dd=ddp) # get the operators
      # project onto the desired low energy manifold
      # now fit the Hamiltonian
      coef = fit_matrix(h,out) # fit the matrix and return dictionary
      try: del coef[("Id")]
      except: pass
      if len(coef)==0: continue # nothing survived for this set; try the next
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
    if isinstance(key,str): return key + "  " # a plain key is not a sequence
    out = ""
    for k in key: out += k + "  "
    return out


def scale2latex(c,tol=1e-12):
    """Format the overall (eV) prefactor of the effective Hamiltonian

    zform is meant for dimensionless ratios: it snaps anything below its own
    absolute tolerance of 1e-3 to "0", which silently erased meV-scale spin
    Hamiltonian parameters. Print the physical scale as a real number instead.
    """
    re,im = float(np.real(c)),float(np.imag(c))
    if np.abs(im)>tol: return "({:.4e}{:+.4e}i)".format(re,im)
    return "{:.4e}".format(re)


def dict2latex(d,tol=1e-4):
    """Transform the dictionary into a latex form"""
    cs = [d[key] for key in d] # coefficients
    cmax = [iy for (ix,iy) in sorted(zip(np.abs(cs),cs))][-1] 
    keys = [key for key in d] # get the keys
    keys = [iy for (ix,iy) in sorted(zip(-np.abs(cs),keys))] # sort the keys
    terms = [] # the terms that survive the tolerance
    for key in keys: # loop
        c = np.round(d[key]/cmax,4) # round the number
        if np.abs(c)<tol: continue
        if .99<np.real(c)<1.01 and np.abs(np.imag(c))<tol: s = "  " # unit coefficient
        else: s = zform(c) + "  " # normalize
        terms.append(s + key2latex(key)) # create the name
    out = "\\begin{aligned}\n"
    # \big instead of \left/\right: a row break (\\) inside a \left...\right
    # group is a hard LaTeX error as soon as three terms survive
    out += "H = \n"+scale2latex(cmax)+" \\big [ " # output string
    for (ik,t) in enumerate(terms): # loop over surviving terms
        out += t
        if ik<len(terms)-1: # separator only between terms
            out += " + \n" # new line
            if (ik+1)%3==0: out += "\\\\ \n" # new line
    out += " \\big ] \n" # last line
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
    else: return False

def braket(a,b):
    return np.abs(np.conjugate(a).dot(b))

from .latexalgebra import matrix2vector



def effective_spin_hamiltonian(lowest,H=None,n=2,nt=2,operators=None):
    """Compute the effective Hamiltonian in Latex form"""
    # get the Hmailtonian
    if H is None: H = lowest.h
    h = lowest.get_representation(H,n=n) # Hamiltonian
    h = h - np.identity(h.shape[0])*np.trace(h)/h.shape[0] # no trace
    atom = lowest.atom # get the atom object
    text = "Hamiltonian written in the low energy manifold with "+str(n)+" states\n"
    if operators is None:
        ops = get_s_operators(atom) # LJ operators
    else:
        ops = dict()
        for key in operators: ops[key] = atom.Operator[key]
    # project to the low energy manifold
    for key in ops:
        O = lowest.get_representation(ops[key],n=n) 
        O = renormalize_spin_operator(O) # renormalize
        ops[key] = O # store
#    exit()
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
    text += fit_warning() # say so if the fit is poor
    return text



def fit_warning():
    """Latex note when the fitted spectrum does not match the real one"""
    err = getattr(fit_matrix,"relative_error",0.)
    if err<=1e-2: return ""
    return ("\n\n\\textbf{Warning:} this effective Hamiltonian reproduces the "
            "spectrum of the manifold only to %.1f\\%%; the operator basis "
            "cannot represent it.\n\n"%(100*err))


def renormalize_spin_operator(m):
    """Rescale a projected spin operator to the pseudo-spin convention

    The operator is scaled so that its largest eigenvalue in magnitude equals
    S = (n-1)/2 for an n-dimensional manifold, i.e. it has the spectrum of a
    spin-S component and the printed \\hat S_x/S_y/S_z labels mean what they
    say. Previously each component was divided by its own smallest |eigenvalue|
    above a hard-coded 0.1, which scaled the transverse and axial components
    differently and made the fitted coefficients incomparable."""
    from .dynamicstk import algebra
    n = m.shape[0] # dimension of the manifold
    S = (n-1)/2. # pseudo-spin
    es = np.abs(algebra.eigvalsh(m)) # magnitudes of the eigenvalues
    emax = np.max(es) if len(es)>0 else 0.
    if emax<1e-8: return m*0. # the operator vanishes in this manifold
    return m*(S/emax)


