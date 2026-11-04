"""Shared conventions: V(x)=s E_R sin^2(k_L x), q=kd, d=pi/k_L.

Static/continuum routines use E_R=1, hbar=1 and tau_R=E_R*t/hbar.
Nearest-neighbor routines instead use J=1 and tau_J=J*t/hbar.
No random numbers or external data are used.
"""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.linalg import eigh_tridiagonal

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "figures"
DATA = ROOT / "results"
TEAL, INK, GRAY, LIGHT = "#1F6F5F", "#171717", "#666666", "#ABABA5"
COLORS = [INK, GRAY, TEAL, "#777777"]
STYLES = ["-", "--", "-", ":"]
CMAP = LinearSegmentedColormap.from_list("probability", ["white", "#DDE9E5", TEAL, "#102C27"])

def setup():
    FIG.mkdir(exist_ok=True, parents=True)
    DATA.mkdir(exist_ok=True, parents=True)
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 9,
        "axes.labelsize": 9, "axes.titlesize": 10, "legend.fontsize": 8,
        "xtick.labelsize": 8, "ytick.labelsize": 8,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": .7, "lines.linewidth": 1.5,
        "legend.frameon": False, "savefig.dpi": 220,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })

def savefig(fig, name):
    """Fixed physical dimensions are set by each caller, not by tight cropping."""
    fig.savefig(FIG / f"{name}.pdf")
    fig.savefig(FIG / f"{name}.png", dpi=220)
    plt.close(fig)

def write_json(name, data):
    (DATA / name).write_text(json.dumps(data, indent=2) + "\n")

def static(q, s, cutoff=12, vectors=False):
    """Bloch Hamiltonian in exp[i(k+2*l*k_L)x]; q=kd=pi*k/k_L.

    The real symmetric tridiagonal matrix is diagonalized without a grid in x.
    cutoff=M keeps l=-M,...,M. Energies are in recoil units.
    """
    ell = np.arange(-cutoff, cutoff + 1)
    diagonal = (q / np.pi + 2 * ell)**2 + s / 2
    off = np.full(2 * cutoff, -s / 4)
    return eigh_tridiagonal(diagonal, off, eigvals_only=not vectors)

def band_fourier(s, cutoff=12, nq=2048, nhop=20):
    """E(q)=Ebar-2 sum_l J_l cos(l*q); all returned energies in E_R."""
    q = np.linspace(-np.pi, np.pi, nq, endpoint=False)
    energy = np.array([static(x, s, cutoff)[0] for x in q])
    mean = float(energy.mean())
    hopping = np.array([-np.mean(energy * np.cos(l*q)) for l in range(1, nhop+1)])
    return q, energy, mean, hopping

def dispersion(q, mean, hopping):
    l = np.arange(1, len(hopping) + 1)
    return mean - 2 * np.cos(np.asarray(q)[..., None] * l) @ hopping
