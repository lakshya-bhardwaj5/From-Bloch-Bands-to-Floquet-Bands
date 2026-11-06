# Numerical methods and provenance

## Source and equation mapping

The source is Chapter 2 of Lakshya Bhardwaj's supplied report, *From Bloch Bands to Floquet Bands: Controlling Atomic Transport in a One-Dimensional Optical Lattice*. The source PDF is identified by its SHA-256 below; the previous repository's report and unrelated calculations have been removed from the current tree.

| Report | Implementation |
|---|---|
| Eqs. (2.3), (2.23)–(2.27) | `lattice.eigensystem`: continuum optical lattice in the minus-reciprocal-vector convention |
| Eqs. (2.28)–(2.29) | Weak-depth edge gap, compared with `s/2` |
| Eq. (2.30) | `lattice.band_metrics`: lowest width and first edge gap |
| Eq. (2.32) | Independent SciPy Mathieu edge checks |
| Eqs. (2.36)–(2.38) | `lattice.fourier_hopping`: Fourier extraction and nearest-neighbor comparison |
| Table 2.1; Figs. 2.1–2.2 | Reproduced numerical table, dispersions, and tight-binding errors |

## Fixed conventions

- Potential: `V(x)/E_R = s sin²(k_L x)`, with `s >= 0` and its mean `s/2` retained.
- Lattice spacing: `d = pi/k_L`; recoil energy: `E_R = hbar² k_L²/(2m)`.
- Momentum: `kappa = k/k_L`; the Fourier coordinate is `q = kd = pi*kappa`.
- Basis: `exp[i(k - 2 ell k_L)x]`, ordered by `ell = -M,...,M`.
- Diagonal entries: `(kappa - 2 ell)² + s/2`; nearest off-diagonal entries: `-s/4`.
- Eigenvector columns have unit sum of squared plane-wave amplitudes. No real-space grid is used to diagonalize the Hamiltonian.
- Band indices start at `n = 0` and follow increasing energy at each momentum. At exactly degenerate free-particle points an individual eigenvector is not uniquely defined; the component animation therefore uses a nonzero lattice depth.

The minus-sign convention follows the current Chapter 2 derivation directly. Reversing `ell` produces the common plus-sign convention without changing eigenvalues.

## Grids and outputs

Static band plots use 401 points from `-1` through `+1` and `M = 12` (25 plane waves). The equal boundary energies are shown twice only for plotting. Hopping extraction uses 2048 half-open, uniformly spaced points in `q` and computes `J_j = -mean[E_0(q) cos(jq)]`. The depth summary uses 121 depths from 0 through 12. The GIF depth sweep uses increments of 0.25; the momentum GIF uses 61 momenta before reversing the sweep.

`data/bands.npz` contains arrays `kappa`, `depths`, and `energies_ER`, with energy shape `(4, 401, 3)`. CSV files expose the four-depth reference table and the 121-depth sweep. PNG and SVG exports share the same plotted values. Both animations use fixed energy axes, a shared indexed color palette, and infinite looping.

The first GIF recomputes stationary spectra while changing depth. The second changes the momentum label at fixed depth and shows the selected energies and normalized components of the lowest eigenstate. Neither GIF represents real-time evolution, Bloch oscillations, or a Floquet calculation.

## Independent checks

`tests/test_lattice.py` uses Python's standard-library `unittest` runner. It checks:

1. The zero-depth spectrum against sorted folded free-particle energies.
2. Three center and three edge energies against independent Mathieu characteristic values. At positive Mathieu parameter `Q=s/4`, center energies are `s/2 + [a0(Q), b2(Q), a2(Q)]`; edge energies are `s/2 + [b1(Q), a1(Q), b3(Q)]`.
3. The lowest three bands at six depths and 21 momenta per depth against `M=16`, and inversion symmetry under `k -> -k`.
4. The weak-depth gap against `s/2` at `s=0.001`.
5. Every rounded energy and hopping entry in Table 2.1.
6. Eigenvector residuals and orthonormality for the component animation.
7. The normalized nearest-neighbor errors for `s=2` and `s=8`.

These are finite-range numerical checks. They do not assert convergence at arbitrary depth or for high bands. No statistical sampling is used, and no measurement error bars are implied. Exact software versions for the generated assets are in `data/environment.json`; quantitative verification is in `data/validation.json`.

## Physical scope

The model is stationary, one dimensional, noninteracting, and single particle. The first three continuum bands are computed before any isolated-band approximation is made. The tight-binding comparison projects onto the lowest band and drops hopping beyond adjacent cells. Wannier functions themselves are not constructed. Periodic driving, trapping, noise, collisions, and dissipation are outside the calculation.

## Source fingerprint

f923126085b23fb8f6c19bcf37e37701c3de1e1515725142c94ca54d34bdabba
