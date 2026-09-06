# routines to compute the gtensor
import numpy as np
import scipy.linalg as lg


def get_gtensor(atom,h0,tol=None):
    """Return the g-tensor of a Kramers doublet ground state

    The Zeeman term of this code is b.(L + 2S), with b = mu_B B in eV, so the
    magnetic moment operator is M_a = (L + 2S)_a.  Inside a two-fold degenerate
    ground state the projected operators can always be written as

        A_a = sum_b g_ab sigma_b / 2

    for some real matrix g and some (basis dependent) Pauli matrices, i.e. the
    doublet behaves as a pseudo-spin 1/2 with H = b . g . S~.  The individual
    g_ab depend on the arbitrary basis chosen inside the doublet, but

        G_ab = sum_c g_ac g_bc = 2 Tr[A_a A_b]

    does not.  This routine returns g = sqrt(G), the symmetric positive
    semidefinite square root: its eigenvalues are the principal values
    g_x, g_y, g_z and its eigenvectors the magnetic axes.

    The result is exact for the doublet, not a finite difference: the previous
    implementation attempted a numerical second derivative of the lowest
    eigenvalue, which is not differentiable at a degeneracy, and it raised
    unconditionally.
    """
    from . import hamiltonians
    if atom is None:
        raise ValueError("get_gtensor needs the Atom object that built H")
    if tol is None: tol = hamiltonians.tol # degeneracy window, in eV
    evals,evecs = hamiltonians.eigenstates(h0) # full spectrum
    de = np.abs(evals - np.min(evals)) # energies from the ground state
    ndeg = int(np.sum(de<tol)) # degeneracy of the ground state
    if ndeg!=2: # a g-tensor is only defined for a pseudo-spin 1/2
        raise ValueError("The g-tensor needs a two-fold degenerate ground "
                "state (a Kramers doublet for odd ne), but this one is %d-fold "
                "degenerate (tol=%g eV). Lower the symmetry, or set "
                "hamiltonians.tol / the tol keyword to adjust the grouping "
                "window."%(ndeg,tol))
    wfs = [v for (e,v) in zip(evals,evecs) if abs(e-np.min(evals))<tol]
    # magnetic moment operators projected on the doublet
    ms = [atom.lx + 2*atom.sx, atom.ly + 2*atom.sy, atom.lz + 2*atom.sz]
    As = []
    for m in ms:
        A = np.array(hamiltonians.get_representation(wfs,m)) # 2x2 block
        A = A - np.identity(2)*np.trace(A)/2. # traceless part (zero by Kramers)
        As.append(A)
    G = np.zeros((3,3)) # gauge invariant G = g g^T
    for i in range(3):
        for j in range(3):
            G[i,j] = 2.*np.trace(As[i]@As[j]).real
    G = (G + G.T)/2. # symmetrize away the rounding
    es,vs = lg.eigh(G) # principal values are the eigenvalues of G
    es = np.where(es>0.,es,0.) # G is positive semidefinite up to rounding
    g = vs@np.diag(np.sqrt(es))@vs.T # symmetric square root
    return g


def get_principal_g(atom,h0,**kwargs):
    """Return the principal g values and the magnetic axes

    Returns (gs,axes), with gs sorted in increasing order and axes[i] the
    direction associated with gs[i]."""
    g = get_gtensor(atom,h0,**kwargs)
    gs,axes = lg.eigh(g) # symmetric, so this is the principal frame
    return gs,axes.T # one axis per row
