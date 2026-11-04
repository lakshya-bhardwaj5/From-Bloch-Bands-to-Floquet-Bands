# The plane-wave basis and the infinite Bloch matrix

At a fixed Bloch wave number $k$, the unknown wave function is built from many
plane waves. A matrix equation determines their coefficients and the allowed
energies. A three-by-three illustration retains three of those components;
the exact matrix is infinite.

This guide uses the notation of the supplied PHY 581 notes. Report Sections
2.2 and 2.3 give the full derivation and its normalization.

## What the symbols label

| Symbol | Meaning |
|---|---|
| $k$ | Bloch wave number chosen in the first Brillouin zone |
| $K=2\ell k_{\mathrm L}$ | A reciprocal lattice wave number |
| $\ell$ | Integer label of that reciprocal vector, extending through all integers |
| $k-K$ | Wave number of one plane-wave component |
| $n$ | Band index: which eigenstate is selected at the fixed $k$ |
| $c^{(n)}_{k-K}$ | Fourier coefficient of that component in the selected eigenstate |
| $U_G$ | Fourier coefficient of the periodic potential |

The expansion is

$$
\psi_{nk}(x)=\sum_K c^{(n)}_{k-K}e^{i(k-K)x}
=e^{ikx}u_{nk}(x),\qquad u_{nk}(x+d)=u_{nk}(x).
$$

Every reciprocal vector has a corresponding component. For example,
$\ell=-1,0,1$ gives wave numbers $k+2k_{\mathrm L}, k, k-2k_{\mathrm L}$.
All share the same translation phase $e^{ikd}$ over one cell.

## From the potential to the matrix

Write

$$
U(x)=\sum_GU_Ge^{iGx}.
$$

Multiplying a basis wave $e^{i(k-K')x}$ by $e^{iGx}$ produces
$e^{i(k-K'+G)x}$. It contributes to the row labeled by $K$ when
$G=K'-K$. Therefore

$$
H_{KK'}(k)=\frac{\hbar^2(k-K)^2}{2m}\delta_{KK'}+U_{K'-K}.
$$

For $U(x)=V_0\sin^2(k_{\mathrm L}x)$,

$$
U_0=V_0/2,\qquad U_{\pm2k_{\mathrm L}}=-V_0/4,
$$

and the other Fourier coefficients vanish. Only columns whose reciprocal
indices differ by zero or one can contribute to a given row. With
$a_{n\ell}=\sqrt L\,c^{(n)}_{k-2\ell k_{\mathrm L}}$ the equation is

$$
D_\ell a_{n\ell}-g a_{n,\ell-1}-g a_{n,\ell+1}
=\mathcal E_{nk}a_{n\ell},
$$

where $D_\ell=\hbar^2(k-2\ell k_{\mathrm L})^2/(2m)+V_0/2$ and $g=V_0/4$.
The normalized coefficients satisfy $\sum_\ell|a_{n\ell}|^2=1$.

## Why three diagonals do not mean three rows

The equation for $a_{n,1}$ includes $a_{n,2}$. The equation for $a_{n,2}$ includes
$a_{n,3}$, and the same pattern continues on both sides. All these coefficients
belong to the same simultaneous eigenvalue problem.

Tridiagonal means that each row has entries only on the main diagonal and
the immediately adjacent diagonals. It does not prescribe the number of rows.

| Cutoff | Retained reciprocal labels | Matrix size |
|---|---|---|
| $M=1$ | $-1,0,1$ | $3\times3$ |
| $M=2$ | $-2,-1,0,1,2$ | $5\times5$ |
| $M=12$ | $-12,\ldots,12$ | $25\times25$ |

Keeping only $-M\leq\ell\leq M$ discards the endpoint couplings to omitted
components. Increase $M$ until the low-energy bands of interest converge.
The numerical checks do this rather than assume the three-component example
is exact.

## Connecting the lecture notation to the code

`common.static()` uses waves with wave number `k+2*ell_code*k_L`.
The lecture basis uses `k-2*ell*k_L`. Setting `ell_code=-ell` reverses the
ordering of the same basis. It permutes rows and columns without changing
eigenvalues. A numerical equivalence check is recorded in
`validation/notation_equivalence_check.json`.

In this guide, capital $K$ is a reciprocal wave number. In the driven model,
the separately defined dimensionless quantity $K=F_0d/(\hbar\omega)$ is the
shaking strength. The code's `q=kd` is dimensionless, whereas the lecture's
general plane-wave wave number $q$ has units of inverse length.

Each plane wave extends over the whole sample. The later Wannier basis uses
orbitals localized near individual lattice wells; its site index has a
different physical meaning from the reciprocal-vector index here.
