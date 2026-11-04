# Q6–D5 literature and repository evidence map

**Purpose.** Trace each mathematical input in the RGS–Q6–Wilson–D5 programme to (i) an external source, (ii) the internal implementation/proof note, and (iii) the boundary between established mathematics and the still-open RGS bridge. A citation supports only the theorem it actually proves; it does not certify a new Q6-specific calculation or physical interpretation.

## 1. External literature mapped to the claims

| Claim | Primary/authoritative source | What it supports | What it does not support |
|---|---|---|---|
| Compact U(1) character expansion | NIST DLMF §10.25.2 (positive power series) and §10.35.1–10.35.2 (generating/Jacobi–Anger identity): https://dlmf.nist.gov/10.25 and https://dlmf.nist.gov/10.35 ; Kogut & Susskind (1975), Phys. Rev. D 11, 395, https://doi.org/10.1103/PhysRevD.11.395 ; Kogut (1979), Rev. Mod. Phys. 51, 659, https://doi.org/10.1103/RevModPhys.51.659 | Fourier–Bessel expansion, compact-link character orthogonality, electric-flux Hamiltonian framework. | The Q6 incidence matrices, rank-124 fibre, or its Q6-specific Haar identity. |
| Wilson transfer matrix | Creutz (1977), Phys. Rev. D 15, 1128, https://doi.org/10.1103/PhysRevD.15.1128 ; Kogut (1979), DOI above | Relation between Wilson Euclidean action and Kogut–Susskind Hamiltonian by transfer matrix. | A continuum limit or a physical 3+1D reduction for this finite six-cube. |
| Modified-Bessel asymptotics | NIST DLMF §§10.40–10.41: https://dlmf.nist.gov/10.40 and https://dlmf.nist.gov/10.41 | Large-argument and uniform large-order asymptotics. | Uniformity for every flux sector without a specified scaling regime. The repository's fixed-integer-flux expansion is not uniform for arbitrary fluxes growing with beta. |
| Cartan domain IV and Bergman kernel | De Fabritiis (1990), *Differential geometry of Cartan domains of type four*, Rend. Mat. Acc. Lincei 1, 131–138: https://www.bdim.eu/item?id=RLIN_1990_9_1_2_131_0 | Type-IV realization, generic norm, Bergman kernel proportional to the generic norm to power -n, metric geometry. | The Q6 five-mode chart or its determinant audit. |
| Generic norm, Hua integral and polynomial | Roos & Yin (2003), *New classes of domains with explicit Bergman kernel*, arXiv:math/0301322: https://arxiv.org/abs/math/0301322 ; Proposition 2.1 and Theorem 2.3 | Generic-norm integral, convergence half-plane Re(s)>-1, Hua polynomial from domain invariants. | A Wilson-to-Hua pushforward; the paper contains no Q6 lattice or RGS model. |
| Euclidean volume of the Lie ball | Rakhmonov & Abdullayev (2019), *On volumes of matrix ball of third type and generalized Lie balls*, Vestnik Udmurtskogo Univ. 29(4), 548–557: https://doi.org/10.20537/vm190406 | Unit-radius type-IV volume pi^n/(2^(n-1)n!); for n=5, pi^5/1920. | A physical radius selection or derivation of alpha. |
| Discrete Hodge theory | Lim (2020), *Hodge Laplacians on Graphs*, SIAM Review 62(3), 685–715: https://doi.org/10.1137/18M1223101 | Incidence matrices, algebraic Hodge decomposition, cohomology and graph Hodge Laplacians. | The specific Q6 product spectrum, which is derived from six interval factors in this repository. |
| A5 character theory | Serre (1977), *Linear Representations of Finite Groups*, Graduate Texts in Mathematics 42, Springer. | Character orthogonality and the irreducibility criterion used to certify the five-dimensional A5 character. | The concrete Q6 action matrices and their computed traces, which are internal finite-group calculations. |
| Orbit–stabilizer microhistory identity | Standard finite-group action theorem; e.g. Dummit & Foote, *Abstract Algebra*, 3rd ed., Wiley (2004). | The exact orbit cardinality and factorial denominator identity for labelled edge-flux events. | The exact Wilson Bessel measure: its pair-creation summation indices are not just the permutation orbit. |
| Local S5 semantics and Gödel's ontological proof | Blackburn, de Rijke & Venema, *Modal Logic*, Cambridge (2001); Kanckos & Lethen (2019), *The Development of Gödel's Ontological Proof*, Review of Symbolic Logic: https://www.cambridge.org/core/journals/review-of-symbolic-logic/article/development-of-godels-ontological-proof/728D16FEF40BA2054055E378C879F369 ; Benzmüller & Scott (2025), https://doi.org/10.1007/s00605-025-02078-x | Equivalence frames give S5 semantics; the cited history and formalizations describe Gödel/Scott's higher-order modal argument. | The Q6 relabelling frame does not supply Gödel's positivity, essence or necessary-existence axioms and is not a proof of the ontological conclusion. |
| Gödel incompleteness (distinct result) | Gödel (1931), *Über formal unentscheidbare Sätze...*, Monatshefte für Mathematik und Physik 38, 173–198, https://doi.org/10.1007/BF01700692. | The incompleteness theorem in its original formal-arithmetic setting. | It is not Gödel's ontological proof and gives no direct derivation of a Wilson, Hua or physical measure law. |
| Smith normal form and integer kernels | FLINT, *Hermite and Smith normal forms*: https://flintlib.org/doc/algorithms/linear-algebra/hermite-smith.html ; Cohen (1996), *A Course in Computational Algebraic Number Theory*, Springer | Unimodular diagonalization, invariant factors, lattice quotient structure and exact integer normal forms. | Correctness of the Q6-specific matrix inputs; those require exact internal checks. |
| Moment determinacy on [0,1] | Hausdorff moment theorem; Diaconis & Freedman (2004), *The Markov moment problem and de Finetti's theorem: Part I*, Math. Z. 247, 183–199: https://doi.org/10.1007/s00209-003-0633-9 | A probability measure on [0,1] is uniquely determined by its full moment sequence. | Convergence of the Q6 Wilson measures. |
| Weak convergence from moments | Fréchet–Shohat theorem; e.g. Sodin, *The classical moment problem*: https://webspace.maths.qmul.ac.uk/a.sodin/teaching/moment/clmp.pdf | Under determinacy and the appropriate tightness/moment hypotheses, moment convergence implies weak convergence. | Angular convergence on D5; radial moments alone do not determine angular measure. |
| Liouville theorem | Cieliebak, *Symplectic Geometry* lecture notes, Hamiltonian systems: https://ammann.app.uni-regensburg.de/lehre/2017w_symplectic/cieliebak_symplectic-geometry-a.pdf | Hamiltonian flows preserve the symplectic volume omega^n/n!. | Preservation of N^-5 dmu_L for arbitrary Hamiltonians; this additionally requires X_H(N)=0. |
| Wyler expression and normalization objection | Robertson (1971), Phys. Rev. Lett. 27, 1545: https://doi.org/10.1103/PhysRevLett.27.1545 ; Gilmore (1972), Phys. Rev. Lett. 28, 462: https://doi.org/10.1103/PhysRevLett.28.462 | Historical expression and radius-normalization debate. | An RGS derivation of alpha. Numerical proximity is not a QED derivation. |

## 2. Exact specialisation to type IV in complex dimension five

For D5 = D^(IV)_5, the Cartan invariants are rank r=2, a=3, b=0, genus g=5 and complex dimension n=5. The generic norm is

\[
N(z)=1-2z^\dagger z+|z^Tz|^2.
\]

The Hua polynomial specialises to

\[
\chi_5(s)=(s+1)_4(s+5/2)
=(s+1)(s+2)(s+3)(s+4)(s+5/2),
\quad \chi_5(0)=60.
\]

With normalized Euclidean volume dnu_E,

\[
\int_{D^{IV}_5}N(z)^s\,d\nu_E(z)
=\frac{\chi_5(0)}{\chi_5(s)},\qquad \Re(s)>-1.
\]

The exact pushforward law of t=N(z) under nu_E is

\[
f_N(t)=\frac{20}{3}(1-\sqrt t)^4(1+4\sqrt t+t),\quad 0<t<1.
\]

Its moments are those of the product of independent Beta(1,1), Beta(2,1), Beta(3,1), Beta(4,1), and Beta(5/2,1) variables. Hausdorff determinacy on [0,1] identifies this law. Thus f_N is the **Euclidean pushforward density of the norm**, not the Bergman density.

The Bergman volume pulls back instead as

\[
\Phi^*(dV_B)=N(\Phi(q,p))^{-5}d\mu_L
\]

on the preimage of D5 under the specific Maxwell-to-D5 chart and normalization in the internal theorem. Its truncated radial normalizer is

\[
Z_5(\epsilon)=\int_\epsilon^1t^{-5}f_N(t)\,dt
\sim\frac{5}{3}\epsilon^{-4}.
\]

The coefficient 5/3 is exact. Hua theory does not determine the physical value of epsilon.

## 3. Internal proof and implementation map

Primary public repository: https://github.com/fabianleonaressi95/rgs-qlg-quaternion-bridge

| Internal result | Canonical files | Evidence status |
|---|---|---|
| Q6 incidence complex and Hodge operators | src/q6_complex.py; docs/Q6_COMBINATORIAL_HODGE_LAPLACIANS.md | Exact finite integer construction. Distinguish the 2-skeleton from the full six-cube when discussing higher Betti numbers. |
| Coexact spectrum and lambda=12 five-mode sector | docs/RGS_Q6_WILSON_D5_MAXWELL_MATCHING.md; tests/test_q6_hodge_spectral_theorem.py | Product-complex derivation plus reproducible finite-matrix audit. |
| A5 irreducible five-mode representation | docs/Q6_A5_TOP_COEXACT_CHARACTER_THEOREM.md; src/rgs_q6_a5_character_theorem.py | Finite-dimensional representation theorem; numerical traces audit the projector argument. |
| U(1) Bessel expansion and Gauss law | docs/RGS_Q6_U1_CHARACTER_ELECTRIC_FLUX.md; docs/RGS_Q6_U1_TRANSFER_OPERATOR.md | Exact Fourier identity and integer Gauss constraint; transfer energy asymptotic has a fixed-flux regime. |
| Integer five-mode quotient and saturated rank-124 kernel | docs/Q6_D5_INTEGER_FIBER_QUOTIENT.md; docs/Q6_D5_SATURATED_KERNEL_AND_TAIL_BOUND.md | Exact rational rank/SNF/HNF calculation; convergence does not mean the full fibre sum has been numerically evaluated. |
| Annihilator-Haar representation | docs/Q6_D5_WILSON_FIVE_TORUS_MARGINAL.md; docs/Q6_GLOBAL_TORUS_STRUCTURE.md; src/rgs_q6_d5_wilson_fiber_torus.py | Exact Fourier/Haar identity. For the audited lattice, all 124 nonzero SNF invariant factors are one, so H_L is connected and isomorphic to T^68. QMC output is an estimate, not a certified interval. |
| Local Maxwell-to-D5 Bergman/Liouville pullback | docs/Q6_D5_BERGMAN_LIOUVILLE_PULLBACK_THEOREM.md; docs/Q6_D5_BERGMAN_LIOUVILLE_MEASURE_LAW.md; src/rgs_q6_d5_liouville_measure.py | Exact ten-dimensional linear-chart volume identity on the preimage of D5; not a global map of the 129-dimensional Wilson quotient. |
| Hua radial density and cutoff integral | docs/Q6_D5_HUA_NORMALIZATION_OBSTRUCTION.md; docs/Q6_D5_HUA_CUTOFF_NORMALIZATION.md; src/rgs_q6_d5_hua_cutoff.py | Exact analytic density, finite-cutoff primitive and divergence coefficient; no physical epsilon. |
| Wyler radius dependence and alpha gate | docs/WYLER_D5_INVARIANT_MEASURE_AND_SCALE_THEOREM.md; docs/WYLER_FINE_STRUCTURE_CONSTANT_RGS_THEOREM.md | Exact scaling obstruction and conditional prediction criteria; no derived alpha. |

## 4. Repository inventory and relevance

The linked GitHub account exposed these repositories during this audit:

- **Directly relevant and public:** [rgs-qlg-quaternion-bridge](https://github.com/fabianleonaressi95/rgs-qlg-quaternion-bridge). It contains the Q6, Wilson, integer-kernel, D5, Hua and Wyler audits listed above.
- **Related research, but not independent support for these theorems:** [naressi-golden-resonance](https://github.com/fabianleonaressi95/naressi-golden-resonance) (private; golden-resonance simulation framework), [naressi-zuliani-law-quantum-morphodynamics](https://github.com/fabianleonaressi95/naressi-zuliani-law-quantum-morphodynamics) (private; toy Monte Carlo), [kilaw](https://github.com/fabianleonaressi95/kilaw) (private; mixed quantum/morphodynamic experiments), and [biverse](https://github.com/fabianleonaressi95/biverse) (private; currently a single README-level Yang–Mills-gap assertion, without a proof artifact in the inspected repository tree).
- **Not direct evidence for the Q6–D5 chain:** [quantum-lindblad](https://github.com/fabianleonaressi95/quantum-lindblad), [cosmodex](https://github.com/fabianleonaressi95/cosmodex), [cosmodex-pwa](https://github.com/fabianleonaressi95/cosmodex-pwa), [q-proof](https://github.com/fabianleonaressi95/q-proof), [ASTRAQ](https://github.com/fabianleonaressi95/ASTRAQ) and [astra-q-core](https://github.com/fabianleonaressi95/astra-q-core) (the latter was empty when inspected).

Private repositories are not cited as external validation. A project README or simulation output is not a peer-reviewed theorem and does not replace a proof. In particular, the biverse README statement is not treated as a proved Yang–Mills mass-gap result.

## 5. Proof obligations no external citation can discharge

1. **Wilson-to-Hua pushforward:** construct the full canonical ten-real-dimensional measure from the five-mode electric quotient and its spatial conjugate sector; prove tightness and weak convergence.
2. **Radial and angular identification:** establish the full pushforward, not only five moments of N. The moment problem identifies the radial law on [0,1], not angular correlations.
3. **Bergman weight:** prove that the limiting measure has Radon–Nikodym density N^-5 relative to the specified Euclidean measure, rather than merely radial law f_N.
4. **Regulator:** derive epsilon from the Wilson transfer spectrum/action and physical scaling variables. A finite spectral threshold lambda_n >= tau is not automatically N >= epsilon.
5. **Continuum/physical normalization:** prove the relevant continuum limit, fix radius and field/charge normalization independently, and state the renormalization scale before claiming a prediction for alpha.

**Conclusion.** The standard ingredients—Bessel character expansion, integer harmonic analysis, Smith normal form, finite Hodge theory, type-IV Bergman geometry, Hua's integral, moment determinacy and Liouville's theorem—have authoritative sources. The new RGS content is their exact assembly on Q6. The Wilson–Hua weak-limit theorem and physical cutoff remain open and are not asserted here.
