# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

`tranci` performs configuration-interaction (CI) calculations for a d-shell (transition-metal atom): 10 spin-orbitals, 1-9 electrons, with Coulomb interaction, crystal fields, spin-orbit coupling and Zeeman terms. Two front-ends share one library: a PyQt5 GUI that emits a LaTeX/PDF summary, and direct Python use from scripts/notebooks.

## Commands

There is **no build step, no test suite and no linter**. The scripts in `examples/*/` are the de facto smoke tests:

```bash
cd examples/octahedral && python main.py              # fastest sanity check (~seconds)
cd examples/effective_hamiltonian && python main.py   # exercises the jax-based fitting path
```

Both must be run **from inside the example directory**: they locate `src/` as `os.getcwd() + "/../../"`, not relative to `__file__`, so `python examples/octahedral/main.py` from the repo root fails on import. The notebooks in `notebooks/` resolve the repo the same way (`os.path.abspath(os.path.join(os.getcwd(), ".."))`), so their kernel must be started with the cwd in `notebooks/`.

The library is **not installed as a package**; every example and notebook prepends the source dir:

```python
import sys; sys.path.append("<repo>/src")
from tranci.atom import get_atom
```

GUI: `bin/tranci` on Linux/Mac, `bin\tranci.bat` on Windows (both just run `python interface_pyqt/main.py`). `install.py` only puts `bin/` on the `PATH` — appending to `~/.bashrc`/`~/.bash_profile` on Linux/Mac, writing `HKCU\Environment\Path` via `winreg` on Windows. The GUI chdirs into a fresh `tempfile.mkdtemp()` folder, so its outputs land there; "Save data" copies `*.tex`, `*.pdf`, `*.OUT`, `*.json` into `./tranci_data` relative to where `tranci` was launched. Every run also writes `parameters.json` (all widget values, keyed by widget `objectName`) next to the results, and File → Save/Load parameters round-trips the same JSON through `qtwrap.get_all_inputs`/`set_all_inputs`, so widget names in `interface.ui` are a file format: renaming one silently orphans that key in old parameter files. Integer fields (`n`, `steps`, `nplot`, `nwf_heff`) are `QSpinBox`es and `qtwrap.get` returns their `.value()`; the rest are `QLineEdit`s with a C-locale double validator. The sweep and operator comboboxes are resolved through the `SWEEP` and `get_operator` tables in `main.py`, so a new combobox entry needs a matching table entry. `interface_pyqt/formulas.py` renders the second-quantized operator of every parameter with matplotlib mathtext (Computer Modern fontset, transparent PNG, device-pixel-ratio 2) into the empty `f_<widget>` QLabels declared in `interface.ui`; adding a parameter means adding both the label in the `.ui` and the entry in `FORMULAS` (a label with no entry stays blank, an entry with no label is skipped). Mathtext is a LaTeX subset: no `\mathbf` inside `\hat`, no `\text`, no environments. All user feedback goes through `qtwrap.status` (status bar) and `show_info`/`show_warning`/`show_error` dialogs, never `print` alone; long callbacks run inside the `busy()` context manager, which disables the window and drives the status-bar progress bar. It needs PyQt5 and `pdflatex` (MiKTeX on Windows) on the system.

`interface_pyqt/interface.py` is **generated** — edit `interface.ui` (Qt Designer) and run `interface_pyqt/convert.sh` (`pyuic5`). Never hand-edit `interface.py`.

`update.py` is just `git add . && git commit -m 'New version' && git push` — that is why every commit is titled "New version".

## Regenerating the operator library (`src/tranci/cilib/`)

The `.op` matrices are precomputed by the C++ code in `src/tranci/cpplib/` and are committed, so this is rarely needed:

```bash
cd src/tranci/cpplib/src && ./compile.sh    # builds main.x; NOTE: hardcoded Eigen path /home/jose/apps/eigen-3.2.5
cd .. && python runall.py                   # loops ne=1..9, writes nelectrons.in, fills cilib/<ne>/
```

`compile.sh` is a Linux/Mac-only bash script; `cilib/` is committed, so a Windows user never regenerates it. `runall.py` and `clean.py` themselves are cross-platform and will pick up `main.exe` if one was built.

`runall.py` starts by deleting `cilib/` with `shutil.rmtree`, so run it only from `src/tranci/cpplib/` and expect to move the result into place. It now aborts (`check=True`) if the C++ generator fails, instead of copying the previous occupation's files.

## Architecture

Data flow: **C++ (`cpplib`) → sparse `.op` files in `cilib/<ne>/` → `CIatom` → user-built Hamiltonian → `Lowest_States`.**

- `atom.get_atom(ne)` is the entry point. It reads `cilib/<ne>/` (all matrices already in the many-body basis for that electron count) and builds `Atom.Operator`, a dict of many-body operators: spin (`sx,sy,sz,s2`), orbital (`lx..l2`), total (`jx..j2`), `Coulomb`/`vc`, spin-orbit `ls`, cartesian crystal-field terms (`x2,y2,z2,x4,y4,z4,x2y2`), and cubic-harmonic occupation projectors (`dz2,dxy,dxz,dyz,dx2y2`).
- Users construct the Hamiltonian themselves as a linear combination, e.g. `H = 4*V + 0.3*CF + 0.1*LS`, then `M = Atom.get_manifolds(H)` → a `hamiltonians.Lowest_States` object exposing `get_gs_degeneracy`, `get_excitations`, `get_gs_projected_eigenvalues(op)`, `get_dynamical_correlator`, `get_correlation_entropy`.
- `Atom.SP_Operator` holds the same operators as **10×10 single-particle** matrices (for `ne>1` it is `get_atom(ne=1).Operator`). `Atom.one2many(m)` lifts a 10×10 matrix into the many-body basis via `edtk/states.py:one2many_basis`. This is the route for custom crystal fields — see `examples/single_occ_VS_many_occ`.
- Orbital index convention (from `orbital.in`): indices 0-4 are m = -2..+2 spin-up, 5-9 are m = -2..+2 spin-down. Rows of `basis.out` are the 10-bit occupation vectors in that order, one per many-body basis state.
- `.op` file format: first line `# SIZE = d`, then `i j real imag` rows for a sparse Hermitian matrix. `read.read_matrix` wraps everything in a bare `except` — a malformed or empty file silently yields a **zero matrix** rather than an error.
- Energies are in **eV** everywhere; `templates.py` converts to GHz/meV when writing `SPECTRUM.OUT`.
- Degenerate manifolds are grouped by `hamiltonians.tol` / `ntol`, which are **module-level globals** mutated by the GUI at runtime. `get_gs_degeneracy(tol=None)` returns an **int** — the number of states within `tol` of the minimum — and reads the module global when `tol` is not given (it used to return a Boltzmann weight at a fixed T=1e-4, which is why some callers still wrap it in `int(np.round(deg))`; that is now a no-op). `Lowest_States.evals` holds the **exact** eigenvalues; `ntol` is a grouping window, not a display precision.
- `hamiltonians.build_hamiltonian(atom, p)` exists for the GUI only: `p` must carry `D, E, U, soc, x2y2, trigonal, j` (with `O`, `z4`, `b` optional and defaulting to zero, plus optional `Uc`); dict and object inputs are read through the same accessor, so both forms accept the same keys. Library code bypasses it and assembles `H` directly.
- The Coulomb prefactor is a **dimensionless multiplier** of a fixed Slater-Condon tensor, not `F0` in eV. For `u*Atom.Operator["Coulomb"]` the Racah parameters are `B = 52.9*u meV`, `C = 210*u meV` (verified against the computed d^2 term spectrum: 5B+2C, 15B, 12B+2C, 22B+7C). The GUI field is labelled `U [multiplier]` for this reason.
- `ls` is the **one-electron** spin-orbit constant and is positive for every filling; Hund's third rule (`J=|L-S|` below half filling, `J=L+S` above) comes out of the many-body diagonalization on its own. Do not flip its sign above half filling.
- `gtensor.get_gtensor(atom, h)` returns the 3x3 g-matrix `sqrt(G)`, `G_ab = 2 Tr[A_a A_b]` with `A_a` the ground-doublet block of `(L+2S)_a`; `get_principal_g` returns the principal values and magnetic axes. It is defined **only for a Kramers doublet** and raises otherwise. Reachable as `Lowest_States.get_gtensor()` / `.get_principal_g()`.
- `doc/tranci_manual.tex` (with its PDF) is the user/physics manual. Keep it in step when operator conventions or GUI labels change; build it with `pdflatex` run twice.
- `effectivehamiltonian.py` imports **jax at module load** (pinned to CPU) and fits the low-energy block with a spin/orbital operator basis, returning LaTeX. It is pulled in by `write.write_all` whenever `n` is not `None`; the GUI passes `n` only when the "Fit effective Hamiltonian" checkbox (off by default) is ticked.
- `write.write_all` assembles and writes `spectrum_ci.tex` only; the GUI is what invokes `pdflatex` (twice, for the table of contents). A library user calling `write_all` directly gets a `.tex` and no PDF. It populates `lowest.gs_manifold`/`lowest.manifolds` itself, so it no longer requires the caller to have run `get_gs_manifold()`/`disentangle_manifolds()` first.

## Cross-platform constraints

The code runs on Linux, Mac and Windows; keep it that way when touching anything that talks to the OS:

- No `os.system` with shell commands (`cp`, `rm -rf`, `mkdir`, `> /dev/null 2>&1`). Use `shutil`, `os.makedirs(..., exist_ok=True)`, `glob`, and `subprocess.run([...], stdout=subprocess.DEVNULL)`.
- No hardcoded `/tmp` or `$HOME` — `tempfile.mkdtemp()` and `os.path.expanduser("~")`.
- Opening the PDF is branched on `sys.platform`: `os.startfile` / `open` / `xdg-open`.
- Forward slashes in path *literals* (`path+"sx.op"`, `mainpath+"/../"`) are fine on Windows and are used throughout `atom.py`/`read.py`; leave them alone.
- `effectivehamiltonian.py` imports jax at module load behind a `try` that re-raises an actionable `ImportError`. jax has official Windows CPU wheels, so this is a missing-dependency guard, not a platform exclusion.

## Legacy / dead code — do not edit

- `src/tranci.py` and `src/interface.xml`: the old GTK front-end (`gi.repository`, `../cilib/` paths that no longer resolve). Superseded by `interface_pyqt/`.
- `src/tranci/states.py`: stale copy of `src/tranci/edtk/states.py`, still using the removed `np.int`. Nothing imports it — `atom.py` uses `edtk.states`.
- `src/tranci/fix_coulomb.py`: Python 2 syntax (`print "..."`), does not even compile. Nothing imports it.
