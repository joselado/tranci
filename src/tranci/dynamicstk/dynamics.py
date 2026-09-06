import numpy as np
import scipy.linalg as lg
from . import algebra
from . import kpm
import scipy.sparse.linalg as slg
import inspect

# scipy renamed cg's tolerance keyword tol -> rtol in 1.12 and removed tol in 1.14
_CG_TOL = "rtol" if "rtol" in inspect.signature(slg.cg).parameters else "tol"


def dynamical_correlator(h,A,B,**kwargs):
#    return dynamical_correlator_ED(h,A,B,**kwargs)
#    return dynamical_correlator_kpm(h,A,B,**kwargs)
    return dynamical_correlator_inv(h,A,B,**kwargs)


def dynamical_correlator_kpm(h,A,B,wf0=None,
        delta=3e-2,es=np.linspace(-1.,1,400)):
    if wf0 is None:
        ees,wfs = algebra.lowest_states(h)
        e0,wf0 = ees[0],wfs[0] # low energy
    else:
        e0 = algebra.lowest_states(h)[0][0]
    # get_moments_vivj gives mu_n = <vj|T_n(H)|vi>, so vj must be A^dagger|0>
    # for this to equal <0|A d(E-H) B|0>, as dynamical_correlator_inv computes
    A = np.conjugate(A.T)
    vi = B@wf0 # first wavefunction
    vj = A@wf0 # second wavefunction
    from scipy.sparse import identity
    m = -identity(h.shape[0])*e0+h # matrix to use
    emax = -algebra.lowest_eigenvalues(-m,n=3)[0] # lowest energy
    scale = np.max([np.abs(e0),np.abs(emax)])*3.0
    n = int(2*scale/delta) # number of polynomials
    (xs,ys) = kpm.dm_vivj_energy(m,vi,vj,scale=scale,
                                npol=n*4,ne=n*10,x=es)
    return xs,np.conjugate(ys)*scale/np.pi # return correlator

def dynamical_correlator_inv(h0,A,B,es=np.linspace(-1,10,600),
        wf0=None,e0=None,
        delta=3e-2,mode="cv"):
    """Calculate a correlation function AB in a frequency window"""
    if wf0 is None:
        ees,wfs = algebra.lowest_states(h0)
        e0,wf0 = ees[0],wfs[0] # low energy
    elif e0 is None: # only diagonalise when the caller has not supplied e0
        e0 = algebra.lowest_states(h0)[0][0]
    ## default method
  #  iden = np.identity(h0.shape[0],dtype=np.complex128) # identity
    from scipy.sparse import identity
    iden = identity(h0.shape[0],dtype=np.complex128) # matrix to use
    bwf0 = B@wf0 # loop invariant
    out = []
    for e in es: # loop over energies
        if mode=="full": # using exact inversion
          g1 = algebra.inv(iden*(e+e0+1j*delta)-h0)
          g2 = algebra.inv(iden*(e+e0-1j*delta)-h0)
          g = 1j*(g1-g2)/2. # normalisation by pi is applied once, at the return
          op = A@g@B # operator
          o = algebra.braket_wAw(wf0,op) # correlator
        elif mode=="cv": # correction vector algorithm
            # The correction-vector matrix is even in delta and its right hand
            # side is odd, so x(-delta) = -x(delta) exactly and the two solves
            # of the original implementation collapse into one:
            #     1j*(o(+delta) - o(-delta))/2 = -<wf0|A|x(+delta)>
            o = -solve_cv(h0,wf0,A,bwf0,e+e0,delta=delta) # conjugate gradient
        else: raise # not recognised
        out.append(o)
    return es,np.array(out)/np.pi # return result





def solve_cv(h0,wf0,si,bwf0,w,delta=0.0):
    """Solve the correction vector with the conjugate gradient method

    Returns <wf0|si|x>, where x solves ((h0-w)^2+delta^2) x = -delta*B|wf0>
    and bwf0 is the precomputed B|wf0>."""
    from scipy.sparse import identity
    iden = identity(h0.shape[0],dtype=np.complex128) # matrix to use
    b = -delta*bwf0 # create the b vector
    hw = h0 - w*iden # shifted Hamiltonian, built once
    A = hw@hw + iden*delta*delta # define A matrix
    x,info = slg.cg(A,b,**{_CG_TOL:1e-10}) # solve the equation
    if info!=0: # never silently return a non-converged solve as a physical result
        raise RuntimeError("Conjugate gradient did not converge in solve_cv "
                           "(scipy cg info=%d, w=%s, delta=%s)"%(info,w,delta))
    o = np.dot(np.conjugate(wf0),si@x) # compute the braket
    return o


def dynamical_correlator_ED(h,a0,b0,delta=2e-2,
        es=np.linspace(-1.0,10.0,600)):
    """Compute a dynamical correlator"""
    emu,vs = algebra.eigh(h)
    U = np.array(vs) # matrix
    Uh = np.conjugate(np.transpose(U)) # Hermitian
    # convention: A and B are used as given, matching dynamical_correlator_inv
    A = Uh@a0@U # get the matrix elements
    B = Uh@b0@U # get the matrix elements
    out = 0.0+es*0.0*1j # initialize
    Ad = np.conjugate(A.T)
    Bd = np.conjugate(B.T)
    out = dynamical_sum(emu,es+1j*delta,A,B,out) # perform the summation
    out -= dynamical_sum(emu,es-1j*delta,Bd,Ad,out) # perform the summation
    return (es,-out.imag/(2*np.pi)) # return correlator


try:
    from numba import jit
except ImportError: # numba is optional; fall back to plain python
    def jit(*args,**kwargs):
        if len(args)==1 and callable(args[0]) and not kwargs: return args[0]
        return lambda f: f

@jit(nopython=True)
def dynamical_sum(es,ws,A,B,out):
    """Return the sum giving the dynamical correlator"""
    out = out*0.0 # initialize
    es = es-np.min(es) # remove minimum
    n = len(es) # number of energies
    for iw in range(len(ws)): # loop over frequencies
        i = 0
        for j in range(n): # loop over energies
            tmp = A[i,j]*B[j,i]
            tmp *= 1./(ws[iw]+es[i] - es[j])
            out[iw] = out[iw] + tmp
    return out # return dynamical correlator

