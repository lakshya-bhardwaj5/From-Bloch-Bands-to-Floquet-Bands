"""Chapter 2: static Bloch bands of V(x) = s E_R sin²(k_L x).

The reciprocal basis is exp[i(k - 2 ell k_L)x], exactly as in Eq. (2.27).
Inputs use kappa = k/k_L; energies and hopping amplitudes use E_R.
"""

import numpy as np
from scipy.linalg import eigh_tridiagonal


def eigensystem(kappa, depth, cutoff=12, bands=3, vectors=False):
    """Solve the finite (2M+1)-dimensional Bloch Hamiltonian.

    `ell` is a reciprocal-vector index, whereas `bands` counts energy-sorted
    eigenstates. The mean potential s/2 is retained in all absolute energies.
    Eigenvector columns contain normalized plane-wave amplitudes a[n, ell].
    """
    if not np.isfinite(kappa) or not np.isfinite(depth) or depth < 0:
        raise ValueError("Use finite kappa and a finite nonnegative depth.")
    if not isinstance(cutoff, int) or cutoff < 1:
        raise ValueError("cutoff must be a positive integer.")
    if not isinstance(bands, int) or not 1 <= bands <= 2 * cutoff + 1:
        raise ValueError("bands must fit inside the plane-wave basis.")
    ell = np.arange(-cutoff, cutoff + 1)
    diagonal = (kappa - 2 * ell) ** 2 + depth / 2
    off_diagonal = np.full(2 * cutoff, -depth / 4)
    return eigh_tridiagonal(
        diagonal, off_diagonal, eigvals_only=not vectors,
        select="i", select_range=(0, bands - 1), lapack_driver="stebz",
    )


def dispersion(kappa, depth, cutoff=12, bands=3):
    """Return energies with shape (number of k points, number of bands)."""
    return np.array([eigensystem(float(k), depth, cutoff, bands) for k in kappa])


def fourier_hopping(depth, cutoff=12, points=2048, terms=4):
    """Eq. (2.37): E0(q) = mean - 2 sum_j J_j cos(j q), q = pi*kappa.

    Use a half-open Brillouin zone for the periodic Fourier quadrature.
    Both endpoints may be included in plotting grids, but never counted twice
    in this discrete average.
    """
    q = np.linspace(-np.pi, np.pi, points, endpoint=False)
    energy = dispersion(q / np.pi, depth, cutoff, bands=1)[:, 0]
    hopping = np.array([
        -np.mean(energy * np.cos(j * q)) for j in range(1, terms + 1)
    ])
    return float(energy.mean()), hopping


def band_metrics(depth, cutoff=12):
    """Lowest-band extrema at center/edge and the first edge gap, Eq. (2.30)."""
    center = eigensystem(0.0, depth, cutoff)
    edge = eigensystem(1.0, depth, cutoff)
    return {
        "depth_s": float(depth),
        "E0_center_ER": float(center[0]),
        "E0_edge_ER": float(edge[0]),
        "bandwidth_ER": float(edge[0] - center[0]),
        "edge_gap_ER": float(edge[1] - edge[0]),
    }
