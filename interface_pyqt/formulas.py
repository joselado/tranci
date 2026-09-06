###############################################
#### Second-quantized formula of each      ####
#### parameter of the interface, rendered  ####
#### with matplotlib mathtext into QLabels ####
###############################################

import io
import matplotlib
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
from PyQt5 import QtWidgets
from PyQt5.QtGui import QPixmap

# c^dagger_{m sigma} creates an electron with l_z = m (m = -2..2) and spin
# sigma in the d shell; V_ijkl is the Slater-Condon tensor of Sec. "Coulomb"
# of the manual, and lower case l,s are single-particle operators
CD = r"c^{\dagger}_{m\sigma}" # creation operator
ONE_ORB = CD+r"c_{m'\sigma}" # one-body, spin diagonal
ONE_FULL = CD+r"c_{m'\sigma'}" # one-body, acting on spin too
SUM_ORB = r"\sum_{mm'\sigma}" # sum for a spin-diagonal one-body operator
SUM_FULL = r"\sum_{mm'\sigma\sigma'}" # sum for a general one-body operator

FORMULAS = { # name of the QLabel : mathtext string
  "f_H": r"$\mathcal{H}=U\hat V_{ee}+\hat V_{CF}"
         r"+\lambda\sum_i\vec l_i\cdot\vec s_i"
         r"+\vec B\cdot(\vec L+2\vec S)+\vec J\cdot\vec S$",
  "f_n": r"$\hat N=\sum_{m\sigma}"+CD+r"c_{m\sigma}$",
  "f_U": r"$U\sum_{ijkl}\sum_{\sigma\sigma'}V_{ijkl}\,"
         r"c^{\dagger}_{i\sigma}c^{\dagger}_{j\sigma'}c_{k\sigma'}c_{l\sigma}$",
  "f_soc": r"$\lambda"+SUM_FULL+r"\langle m\sigma|\vec l\cdot\vec s|m'\sigma'\rangle\,"+ONE_FULL+"$",
  "f_D": r"$D\sum_{m\sigma}m^{2}\,"+CD+r"c_{m\sigma}$",
  "f_E": r"$E"+SUM_ORB+r"\langle m|l_x^{2}-l_y^{2}|m'\rangle\,"+ONE_ORB+"$",
  "f_O": r"$O"+SUM_ORB+r"\langle m|l_x^{4}+l_y^{4}+l_z^{4}|m'\rangle\,"+ONE_ORB+"$",
  "f_trigonal": r"$t"+SUM_ORB+r"\langle m|(\vec l\cdot\hat n)^{2}|m'\rangle\,"+ONE_ORB+
                r",\;\;\hat n=(1,1,1)/\sqrt{3}$",
  "f_z4": r"$D'\sum_{m\sigma}m^{4}\,"+CD+r"c_{m\sigma}$",
  "f_x2y2": r"$E'"+SUM_ORB+r"\langle m|(l_xl_y)^{2}+(l_yl_x)^{2}|m'\rangle\,"+ONE_ORB+"$",
  "f_B": r"$\vec B\cdot"+SUM_FULL+r"\langle m\sigma|\vec l+2\vec s|m'\sigma'\rangle\,"+ONE_FULL+"$",
  "f_dir_b": r"$\hat B=(\sin\pi\theta_B\cos\pi\phi_B,\;\sin\pi\theta_B\sin\pi\phi_B,\;\cos\pi\theta_B)$",
  "f_j": r"$\vec J\cdot\sum_{m\sigma\sigma'}\vec s_{\sigma\sigma'}\,"+CD+r"c_{m\sigma'}"+
         r",\;\;\vec s=\frac{\vec\sigma}{2}$",
  "f_dir_j": r"$\hat J=(\sin\pi\theta_J\cos\pi\phi_J,\;\sin\pi\theta_J\sin\pi\phi_J,\;\cos\pi\theta_J)$",
}


def render(tex,color="black",size=11.,scale=2):
  """Render a mathtext string into a QPixmap with transparent background

  The image is rendered at `scale` times the screen resolution and marked
  with that device pixel ratio, so it is crisp on HiDPI screens too."""
  dpi = int(96*scale)
  with matplotlib.rc_context({"mathtext.fontset": "cm"}):
    fig = Figure(figsize=(0.01,0.01),dpi=dpi)
    FigureCanvasAgg(fig) # attach a renderer that does not need a display
    fig.text(0,0,tex,fontsize=size,color=color)
    buf = io.BytesIO()
    fig.savefig(buf,format="png",dpi=dpi,transparent=True,
                bbox_inches="tight",pad_inches=0.02)
  pix = QPixmap()
  pix.loadFromData(buf.getvalue(),"PNG")
  pix.setDevicePixelRatio(scale)
  return pix


def add_formulas(form,size=11.):
  """Fill every f_* QLabel of the interface with its rendered formula.
  If rendering fails the raw TeX string is shown instead."""
  color = form.palette().windowText().color().name() # follow the theme
  for name,tex in FORMULAS.items():
    label = getattr(form,name,None)
    if label is None: continue # not in this interface
    label.setToolTip(tex.strip("$")) # the TeX source, for copy-paste
    try: label.setPixmap(render(tex,color=color,size=size))
    except Exception as e: # mathtext could not parse it
      print("Could not render",name,":",e)
      label.setText(tex)
