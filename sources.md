# Source verification and claim crosswalk

Checked on 2026-10-31 (America/Phoenix). The project does not claim novelty.
No full papers or experimental datasets are redistributed.

## Exact principal texts checked

1. Martin Holthaus, *Floquet engineering with quasienergy bands of periodically
   driven optical lattices*, J. Phys. B 49, 013001 (2016),
   DOI: [10.1088/0953-4075/49/1/013001](https://doi.org/10.1088/0953-4075/49/1/013001).
   Primary full text checked: [arXiv:1510.09042v1](https://arxiv.org/pdf/1510.09042v1),
   submitted 30 October 2015; manuscript dated 26 August 2015. The paper's
   published 2016 year is used in the report citation. Relevant portions:
   Sec. II (cosine-lattice bands), Sec. III (Floquet formalism), Sec. IV B
   (shaking transformations and propagator calculation), Sec. IV C (interband
   resonances). Its Table I gives W0/E_R=0.345, gap=1.969 for s=4, and
   W0/E_R=0.123, gap=3.770 for s=8. Our independently computed values agree
   at that published rounding. This comparison authenticates the conventions
   and finite static calculation, not every Floquet approximation.

2. H. Lignier et al., *Dynamical control of matter-wave tunneling in periodic
   potentials*, Phys. Rev. Lett. 99, 220403 (2007),
   DOI: [10.1103/PhysRevLett.99.220403](https://doi.org/10.1103/PhysRevLett.99.220403).
   Primary full text checked: [arXiv:0707.0403v1](https://arxiv.org/pdf/0707.0403v1),
   3 July 2007. Eq. (2): Bessel hopping renormalization. Fig. 2 and associated
   text: expansion-based |J_eff/J| and suppression of at least a factor 25
   within resolution. Fig. 3 and associated text: time-of-flight coherence,
   shifted interference pattern in negative-hopping region, dephasing and
   excited-band limitations. We do not digitize the experimental points or
   infer exact zero hopping from the measured suppression.

3. A. Eckardt et al., *Exploring dynamic localization with a Bose-Einstein
   condensate*, Phys. Rev. A 79, 013611 (2009),
   DOI: [10.1103/PhysRevA.79.013611](https://doi.org/10.1103/PhysRevA.79.013611).
   Primary full text checked: [arXiv:0812.1997v2](https://arxiv.org/pdf/0812.1997v2),
   11 March 2009. Eq. (18): isolated-band period average; Eq. (22): group
   velocity; Sec. III and Table I: longer-range hopping and waveform dependence;
   Eqs. (30), (33)-(34): Bessel/sinc factors; Sec. IV, Figs. 3-4: expansion
   and dephasing measurements. Their hopping coefficients c_l are -J_l/E_R
   in our convention. Their Gaussian-width parameter is not silently equated
   with our probability standard deviation: our width formula is derived with
   the explicitly normalized variance convention in the report.

## Source-dependent implication and notation checks

| Established result | Required hypotheses | Translation and verdict |
|---|---|---|
| J_eff=J J0(K) | Homogeneous noninteracting nearest-neighbor chain; sinusoidal force | Verified after translating force-energy parameter K in Lignier to F0*d; their K0 is our dimensionless K |
| epsilon=period average of E(k+A/hbar) | Isolated band; negligible interband transfer; appropriate phase/gauge choice | Verified with retraced 1D sinusoidal path; conditional as a claim about a full continuum lattice |
| vbar=(1/hbar) d epsilon/dk | Smooth normalized Floquet branch; narrow packet or weighted average | Derived using extended-space Hellmann-Feynman and stated packet assumptions |
| Flat band implies recurrence | Wave packet supported in that Floquet band; periodic micromotion | Verified; exact in the NN model; not asserted for a generic cosine potential |
| Square forcing has all-harmonic collapse at K=2r | Isolated band, uniform square-wave force, no interband/interaction effects | Verified by explicit impulse integral; no square-wave numerical execution claimed |

Our s*sin^2 convention is Holthaus's +(V0/2)*cos(2*kL*x) after a half-cell
translation and an additive V0/2. Holthaus's beta=F0/(hbar*omega*kL) is K/pi.
We retain the A^2/(2m) scalar term in continuum calculations. Laboratory-frame
quasienergies differ by the stated mean quiver-energy offset; slopes do not.

## Supporting primary publications

- J. H. Shirley, Phys. Rev. 138, B979 (1965),
  [10.1103/PhysRev.138.B979](https://doi.org/10.1103/PhysRev.138.B979).
  Publication metadata and primary full-text availability verified; historical
  attribution to the quantum time-periodic eigenproblem.
- H. Sambe, Phys. Rev. A 7, 2203 (1973),
  [10.1103/PhysRevA.7.2203](https://doi.org/10.1103/PhysRevA.7.2203).
  Primary publication metadata/abstract verified; extended-space quasienergy
  formalism. The report derives its matrix rather than asserting an unexamined theorem.
- D. H. Dunlap and V. M. Kenkre, Phys. Rev. B 34, 3625 (1986),
  [10.1103/PhysRevB.34.3625](https://doi.org/10.1103/PhysRevB.34.3625).
  Primary metadata plus detailed attribution in Eckardt 2009 checked; the
  exact NN localization result is independently derived in the report.
- M. Holthaus, Phys. Rev. Lett. 69, 351 (1992),
  [10.1103/PhysRevLett.69.351](https://doi.org/10.1103/PhysRevLett.69.351).
  Citation/title verified in Holthaus 2016 reference [79]. Historical
  attribution; no new source-dependent technical claim beyond our derivation.
- K. W. Madison et al., Phys. Rev. Lett. 81, 5093 (1998),
  [10.1103/PhysRevLett.81.5093](https://doi.org/10.1103/PhysRevLett.81.5093).
  Citation/title in Holthaus 2016 reference [84], and historical observation
  discussed in Eckardt 2009 introduction. No numerical experimental comparison.
- D. Jaksch et al., Phys. Rev. Lett. 81, 3108 (1998),
  [10.1103/PhysRevLett.81.3108](https://doi.org/10.1103/PhysRevLett.81.3108).
  Primary abstract/full text verified for optical-lattice Bose-Hubbard modeling.
- I. Bloch, J. Dalibard and W. Zwerger, Rev. Mod. Phys. 80, 885 (2008),
  [10.1103/RevModPhys.80.885](https://doi.org/10.1103/RevModPhys.80.885).
  Primary review metadata checked for general optical-lattice quantum simulation.
- S. Arlinghaus and M. Holthaus, Phys. Rev. B 84, 054301 (2011),
  [10.1103/PhysRevB.84.054301](https://doi.org/10.1103/PhysRevB.84.054301).
  Citation/title checked in Holthaus 2016 reference [71]; the exact velocity
  relation needed here is derived directly and also checked in Eckardt 2009.

## Independent numerical additions and non-claims

All figures are newly executed calculations of established models. Specific
parameter comparisons, convergence errors and the two-frequency resonant
excitation example are educational computations, not novel physics claims.
There are no experimental fits, digitized data, simulated time-of-flight
interference images, many-body calculations or indefinite-time stability claims.
