"""Generate the figures used in the README.

Run from this directory:  python make_figures.py
"""
import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from tranci.atom import get_atom

plt.rcParams.update({"font.size": 12, "axes.spines.top": False,
                     "axes.spines.right": False})
here = os.path.dirname(os.path.abspath(__file__))
colors = ["#1f77b4", "#d62728", "#2ca02c", "#ff7f0e", "#9467bd"]

# ---------------------------------------------------------------------------
# 1. Multiplet spectrum of a d^2 ion versus octahedral crystal field
# ---------------------------------------------------------------------------
Atom = get_atom(ne=2)
V = Atom.Operator["Coulomb"]
Oct = Atom.Operator["x4"] + Atom.Operator["y4"] + Atom.Operator["z4"]
u = 2.0                                  # B = 106 meV, C = 420 meV
Os = np.linspace(0.0, 0.5, 120)          # 10Dq = 6 O, in eV
evals = []
for O in Os:
    M = Atom.get_manifolds(u*V + O*Oct)
    evals.append(M.evals)
evals = np.array(evals)
fig, ax = plt.subplots(figsize=(6, 4.2))
ax.plot(6*Os, evals, color="#1f77b4", lw=1.0, alpha=0.8)
ax.set_xlabel(r"Octahedral splitting $10Dq$ [eV]")
ax.set_ylabel(r"$E - E_0$ [eV]")
ax.set_ylim(-0.05, 4.5)
ax.set_xlim(0, 6*Os[-1])
ax.set_title(r"$d^2$ multiplets: Coulomb + octahedral field", fontsize=12)
ax.text(0.02, 0.05, r"$^3F$ (21-fold)", transform=ax.transAxes, color="#1f77b4")
fig.tight_layout()
fig.savefig(os.path.join(here, "d2_multiplets_vs_10Dq.png"), dpi=160)

# ---------------------------------------------------------------------------
# 2. Magnetic anisotropy: Zeeman splitting of a d^7 ion in a uniaxial field
# ---------------------------------------------------------------------------
Atom = get_atom(ne=7)
V = Atom.Operator["Coulomb"]
D = Atom.Operator["z2"]
LS = Atom.Operator["ls"]
Lz, Sz = Atom.Operator["lz"], Atom.Operator["sz"]
Lx, Sx = Atom.Operator["lx"], Atom.Operator["sx"]
H0 = 2.0*V + 0.1*D + 0.05*LS
muB = 5.7884e-5                          # eV/T
Bs = np.linspace(0, 20, 60)
ez, ex = [], []
for B in Bs:
    b = muB*B
    ez.append(Atom.get_manifolds(H0 + b*(Lz + 2*Sz)).evals[:6])
    ex.append(Atom.get_manifolds(H0 + b*(Lx + 2*Sx)).evals[:6])
ez, ex = 1e3*np.array(ez), 1e3*np.array(ex)
fig, axs = plt.subplots(1, 2, figsize=(8, 3.8), sharey=True)
for ax, e, lab in zip(axs, [ez, ex], [r"$\vec B \parallel z$", r"$\vec B \parallel x$"]):
    ax.plot(Bs, e, color="#d62728", lw=1.5)
    ax.set_xlabel("Magnetic field [T]")
    ax.set_title(lab, fontsize=12)
    ax.set_xlim(0, Bs[-1])
axs[0].set_ylabel(r"$E - E_0$ [meV]")
fig.suptitle(r"$d^7$ ion, uniaxial crystal field + spin-orbit: magnetic anisotropy", fontsize=12)
fig.tight_layout()
fig.savefig(os.path.join(here, "d7_zeeman_anisotropy.png"), dpi=160)

# ---------------------------------------------------------------------------
# 3. Spin crossover of a d^6 ion (Fe2+) versus octahedral crystal field
# ---------------------------------------------------------------------------
Atom = get_atom(ne=6)
V = Atom.Operator["Coulomb"]
Oct = Atom.Operator["x4"] + Atom.Operator["y4"] + Atom.Operator["z4"]
t2g = Atom.Operator["dxy"] + Atom.Operator["dxz"] + Atom.Operator["dyz"]
eg = Atom.Operator["dz2"] + Atom.Operator["dx2y2"]
Os = np.linspace(0.0, 0.6, 121)
nt, ne_, S, deg = [], [], [], []
for O in Os:
    M = Atom.get_manifolds(2.0*V + O*Oct)
    nt.append(np.mean(M.get_gs_projected_eigenvalues(t2g)))
    ne_.append(np.mean(M.get_gs_projected_eigenvalues(eg)))
    s2 = np.mean(M.get_gs_projected_eigenvalues(Atom.Operator["s2"]))
    S.append((-1 + np.sqrt(1 + 4*s2))/2)
    deg.append(M.get_gs_multiplicity())
fig, axs = plt.subplots(1, 2, figsize=(8, 3.8))
axs[0].plot(6*Os, nt, color="#1f77b4", lw=2, label=r"$t_{2g}$")
axs[0].plot(6*Os, ne_, color="#ff7f0e", lw=2, label=r"$e_g$")
axs[0].set_xlabel(r"Octahedral splitting $10Dq$ [eV]")
axs[0].set_ylabel("Ground-state occupation")
axs[0].set_ylim(-0.2, 6.3)
axs[0].legend(frameon=False)
axs[1].plot(6*Os, S, color="#d62728", lw=2)
axs[1].set_xlabel(r"Octahedral splitting $10Dq$ [eV]")
axs[1].set_ylabel("Total spin $S$")
axs[1].set_ylim(-0.1, 2.3)
axs[1].text(0.05, 0.85, "high spin\n$^5T_{2g}$", transform=axs[1].transAxes)
axs[1].text(0.7, 0.12, "low spin\n$^1A_{1g}$", transform=axs[1].transAxes)
for ax in axs: ax.set_xlim(0, 6*Os[-1])
fig.suptitle(r"$d^6$ ion: high-spin to low-spin crossover", fontsize=12)
fig.tight_layout()
fig.savefig(os.path.join(here, "d6_spin_crossover.png"), dpi=160)
print("d6 crossover at 10Dq =", 6*Os[np.argmax(np.array(S) < 1)], "eV; degeneracies", sorted(set(deg)))
print("done")
