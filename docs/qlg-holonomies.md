# QLG Holonomies: Quaternionic Lattice Connections

## 1. Purpose

This document describes the proposed connection between quaternionic rotations, SU(2) lattice variables and the QLG discrete geometric framework developed within the broader RGS–PSU research programme.

The immediate objectives are:

1. Define SU(2) link variables on an oriented discrete graph.
2. Construct parallel transport along lattice paths.
3. Define closed-loop and plaquette holonomies.
4. Establish the correct inverse and orientation conventions.
5. Distinguish SU(2) information from its SO(3) projection.
6. Specify reproducible numerical validation procedures.

The mathematical construction is separated from its possible physical interpretation.

## 2. Discrete graph and link variables

Let

\[
G=(V,E)
\]

be an oriented graph, with vertex set \(V\) and oriented edge set \(E\).

Associate to every oriented edge \(g\rightarrow h\) a group element

\[
U_{gh}\in SU(2).
\]

Using the quaternionic representation,

\[
U_{gh}=U(q_{gh}),
\qquad
q_{gh}\in S^3.
\]

The reverse orientation is defined by

\[
U_{hg}=U_{gh}^{-1}=U_{gh}^{\dagger}.
\]

For unit quaternions,

\[
q_{hg}=q_{gh}^{*}.
\]

This relation is mandatory: the reverse link is the group inverse, not an arbitrary antipodal representative.

## 3. Parallel transport

Consider an oriented path

\[
\gamma=(v_0,v_1,\ldots,v_n).
\]

Define the parallel transport along this path by

\[
T_\gamma=
U_{v_0v_1}U_{v_1v_2}\cdots U_{v_{n-1}v_n}.
\]

The product order is fixed by the repository convention.

For a sequence of active transformations, the rightmost factor acts first.

For a closed path,

\[
v_n=v_0,
\]

the parallel transport becomes the holonomy

\[
H_\gamma=
U_{v_0v_1}U_{v_1v_2}\cdots U_{v_{n-1}v_0}.
\]

The holonomy is an element of SU(2).

Its SO(3) projection is

\[
R_\gamma=\rho(H_\gamma).
\]

The two objects contain related but different information.

## 4. Plaquette holonomy

For an oriented plaquette with vertices

\[
g\rightarrow h\rightarrow k\rightarrow l\rightarrow g,
\]

define

\[
H_{ghkl}=
U_{gh}U_{hk}U_{kl}U_{lg}.
\]

The plaquette holonomy measures the net group transport around the closed loop.

For a flat connection on a simply connected region, the holonomy around contractible loops is the identity, subject to the usual assumptions about the discrete connection.

In general,

\[
H_{ghkl}\ne I.
\]

Its conjugacy class is invariant under local gauge transformations.

For SU(2), a useful gauge-invariant quantity is

\[
W_{ghkl}=\frac12\operatorname{Tr}(H_{ghkl}).
\]

For an SU(2) matrix,

\[
W_{ghkl}\in[-1,1].
\]

If the eigenvalues of the holonomy are

\[
e^{i\alpha},\quad e^{-i\alpha},
\]

then

\[
W_{ghkl}=\cos\alpha.
\]

The angle \(\alpha\) is the SU(2) eigenangle. It must not be confused automatically with the full SO(3) rotation angle, which is defined modulo the double-cover identification.

## 5. Gauge transformations

Let each vertex \(g\) carry a local gauge transformation

\[
G_g\in SU(2).
\]

The link variable transforms as

\[
U_{gh}\longmapsto
U'_{gh}=G_gU_{gh}G_h^{-1}.
\]

For a closed loop based at \(g\), the holonomy transforms as

\[
H_\gamma\longmapsto
H'_\gamma=G_gH_\gamma G_g^{-1}.
\]

Therefore,

\[
\operatorname{Tr}(H'_\gamma)
=
\operatorname{Tr}(H_\gamma).
\]

The trace of a closed-loop holonomy is gauge invariant.

This property provides one of the most important independent tests of the lattice implementation.

## 6. Orientation reversal

For a path

\[
\gamma=(v_0,v_1,\ldots,v_n),
\]

the reversed path is

\[
\gamma^{-1}=(v_n,v_{n-1},\ldots,v_0).
\]

Its holonomy is

\[
H_{\gamma^{-1}}=H_\gamma^{-1}.
\]

For unitary SU(2) matrices,

\[
H_\gamma^{-1}=H_\gamma^\dagger.
\]

For unit quaternion representatives,

\[
q_{\gamma^{-1}}=q_\gamma^*.
\]

The antipodal element

\[
-q_\gamma^*
\]

has the same SO(3) image as the inverse, but is not the inverse SU(2) element.

This distinction must be preserved when calculating products around oriented plaquettes.

## 7. SU(2) lifts and SO(3) projection

The projection

\[
\rho:SU(2)\rightarrow SO(3)
\]

satisfies

\[
\rho(U)=\rho(-U).
\]

Consequently, an SO(3) lattice calculation cannot, by itself, reconstruct the individual SU(2) lift signs.

For a closed path, replacing link representatives by

\[
U_j\longmapsto s_jU_j,
\qquad s_j\in\{+1,-1\},
\]

changes the SU(2) holonomy by

\[
H_\gamma\longmapsto
\left(\prod_j s_j\right)H_\gamma.
\]

Its SO(3) projection remains unchanged.

The SU(2) trace, however, changes sign when the total sign product is negative.

Therefore, the implementation must retain the SU(2) link data when lift information is relevant.

## 8. Two-sheet construction

The broader QLG research programme includes a proposed two-sheet or biverse construction.

One computational convention previously considered is to define the second sheet through complex conjugation:

\[
U^{(-)}_{gh}=\overline{U^{(+)}_{gh}}.
\]

This is a specific construction and must not be confused with the quaternionic inverse or with an independently chosen antipodal lift.

Complex conjugation, Hermitian conjugation and quaternionic conjugation are different operations:

\[
\overline U,\qquad U^\dagger,\qquad q^*.
\]

Their relationships must be derived using the chosen matrix and quaternion representations.

For the present implementation, the two sheets are treated as separate datasets until their transformation properties and physical interpretation have been independently established.

No physical parity, time-reversal or particle interpretation follows from complex conjugation alone.

## 9. Previously reported QLG computational structure

The following quantities have been reported in previous QLG computational work and are included here as reference targets for future independent reproduction.

| Quantity | Previously reported value |
|---|---:|
| Edge rotations | 192 |
| SU(2) lifts | 192 |
| Plaquette holonomies | 240 |
| SU(2) unitarity residual | \(3.21\times10^{-16}\) |
| SU(2) to SO(3) lift residual | \(9.49\times10^{-16}\) |
| SO(3) orthogonality residual | \(1.17\times10^{-15}\) |
| Wilson angle RMS | 1.82391 |
| Wilson angle maximum | 5.90535 |

These values are historical numerical reports from the broader project.

They are not results independently reproduced by this repository and must not be labelled as verified here until the relevant input data, algorithms and tests are available.

In particular, the definition of the reported Wilson angle must be documented before comparing it with the SU(2) eigenangle or the SO(3) rotation angle.

## 10. Numerical validation protocol

The initial QLG implementation should include the following tests.

### 10.1 Link validity

For every link,

\[
U_{gh}^{\dagger}U_{gh}=I,
\]

\[
\det U_{gh}=1.
\]

### 10.2 Reverse-link consistency

Verify

\[
U_{hg}U_{gh}=I.
\]

### 10.3 Holonomy validity

For every closed path,

\[
H_\gamma^\dagger H_\gamma=I,
\]

\[
\det H_\gamma=1.
\]

### 10.4 Orientation reversal

Verify

\[
H_{\gamma^{-1}}=H_\gamma^{-1}.
\]

### 10.5 Gauge covariance

Generate random vertex transformations and verify

\[
H'_\gamma=G_{v_0}H_\gamma G_{v_0}^{-1}.
\]

### 10.6 Gauge-invariant trace

Verify

\[
\operatorname{Tr}(H'_\gamma)
=
\operatorname{Tr}(H_\gamma).
\]

### 10.7 SU(2)/SO(3) consistency

Verify

\[
\rho(H_\gamma)
=
\rho(U_1)\rho(U_2)\cdots\rho(U_n).
\]

### 10.8 Lift-sign consistency

Explicitly test the effect of changing one or more link signs and verify that the SU(2) product changes as expected while its SO(3) projection remains invariant.

Every test must report its residual, input dimensions and numerical tolerance.

## 11. Open mathematical questions

The following questions remain open within the broader QLG construction:

1. What conditions guarantee a globally consistent choice of SU(2) lifts on the discrete graph?
2. How should central sign information be organized around independent cycles?
3. What is the precise relationship between the two-sheet construction and lattice gauge transformations?
4. Can the discrete connection be related to a continuum connection through a controlled convergence procedure?
5. Under what assumptions can the lattice holonomies converge to continuum Wilson loops?
6. What additional structure is required to derive a physically meaningful gauge theory?

These questions require explicit definitions and independent proofs.

## 12. Relationship to RGS–PSU

The QLG construction is intended as a mathematical bridge between discrete geometric structures and group-valued parallel transport.

The broader RGS–PSU programme also investigates:

- discrete geometric operators;
- tangent-space reductions;
- quaternionic and spinorial representations;
- teleparallel structures;
- spectral convergence;
- possible continuum gauge-theoretic limits.

These connections must be established through explicit maps and convergence arguments.

In particular, a verified SU(2) lattice connection does not automatically establish a Yang–Mills continuum limit or a physical identification with a Standard Model gauge sector.

## 13. Scientific status

**Established mathematical identities**

- SU(2) link inversion.
- Closed-loop holonomy composition.
- Gauge covariance of holonomies.
- Gauge invariance of the trace.
- Orientation reversal by group inversion.
- The SU(2)/SO(3) double-cover relation.

**Previously reported numerical evidence**

- The 192-edge and 240-plaquette QLG calculations.
- The residuals listed in Section 9.

**Open research questions**

- Global lift consistency.
- Two-sheet interpretation.
- Continuum convergence.
- Physical gauge-theory identification.

The purpose of this document is to provide a precise foundation for independent implementation and review.
