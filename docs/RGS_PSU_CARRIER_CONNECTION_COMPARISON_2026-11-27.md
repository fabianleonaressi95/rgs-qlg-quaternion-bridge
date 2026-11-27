# RGS–PSU archived connection versus carrier-derived coframe

**Date:** 27 November 2026  
**Classification:** exact generator-level connection equivalence; finite-Q6 coframe and torsion comparison implemented as executable tests. The original A13 NPZ binary payload is not loaded by this audit.

## 1. The connection is the same construction in two frames

The A13 builder defines, on each positive Q6 link \(g\to h\),
\[
n_g=\frac{X_g}{\|X_g\|},\qquad
R^{\rm A13}_{gh}=\operatorname{Rot}_{\min}(n_g,n_h),
\]
and lifts it to \(U^{\rm A13}_{gh}\in SU(2)\) with the scalar-first quaternion branch. The carrier tetrahedron audit defines the same ambient shortest rotation \(R^g_{gh}\), then expresses it in deterministic radial frames \(F_g,F_h\in SO(3)\):
\[
R^{\rm car}_{gh}=F_h^T R^g_{gh}F_g.
\]
Since \(R^g_{gh}=R^{\rm A13}_{gh}\) by the two generator definitions,
\[
\boxed{R^{\rm car}_{gh}=F_h^T R^{\rm A13}_{gh}F_g.}
\]
This is not a numerical fit: it is the exact local-frame gauge transformation of the same edge connection.

If \(S_g,S_h\in SU(2)\) lift \(F_g,F_h\), then the corresponding spin link is
\[
U^{\rm car}_{gh}=S_h^\dagger U^{\rm A13}_{gh}S_g,
\]
up to the unavoidable \(\pm1\) choice in the double cover. The audit checks the adjoint \(SU(2)\to SO(3)\) image rather than comparing lift signs.

The edge ordering is the A13/RGS vertex-major, least-significant-bit ordering. The already merged Q6 alignment permutation maps it bijectively to the canonical lexicographic/big-endian edge ordering; the edge-set identity and SHA-256 are checked again by the executable audit.

## 2. Structural reduction in radial frames

Because \(R^{\rm A13}_{gh}n_g=n_h\) and the third column of \(F_g\) is \(n_g\),
\[
R^{\rm car}_{gh}(0,0,1)^T=(0,0,1)^T.
\]
Therefore each local link lies in the \(SO(2)\) stabilizer of the radial axis. Its spin lift lies in the corresponding diagonal \(U(1)\subset SU(2)\), up to sign. This is a property of the specified radial-frame gauge, not a claim that the full PSU theory is an Abelian gauge theory.

## 3. The coframes are not identical

The archived A13 feature is the normalized four-component edge vector
\[
\widehat e^{\rm A13}_{gh}
=\frac{(\Delta X_{gh},\Delta r_{gh})}
{\sqrt{\|\Delta X_{gh}\|^2+(\Delta r_{gh})^2}},
\qquad \Delta r_{gh}=\|X_h\|-\|X_g\|.
\]
Its fourth component is a radial scalar. It is not Lorentz time.

The carrier construction instead uses a three-component spatial coframe in the local radial frame:
\[
E^{\rm car}_{gh}
=\frac12\left(F_g^T\Delta X_{gh}
+(R^{\rm car}_{gh})^T F_h^T\Delta X_{gh}\right).
\]
It is a symmetric, gauge-covariant average of the source chord and the transported target chord. It is not obtained by merely dropping the fourth A13 component: the normalization and symmetrization are different. The script reports the edgewise norm ratios and the difference from the archived normalized spatial projection.

The reverse-edge law
\[
E^{\rm car}_{hg}=-R^{\rm car}_{gh}E^{\rm car}_{gh}
\]
and the bound \(\|E^{\rm car}_{gh}\|\leq\|\Delta X_{gh}\|\) are checked on all 192 Q6 links.

## 4. Square torsion comparison

The audit reconstructs the 240 canonical Q6 squares. It computes the covariant face sum twice:
- with the raw local spatial coordinate cochain \(F_g^T\Delta X_{gh}\);
- with the new symmetrized carrier coframe \(E^{\rm car}_{gh}\).

The first is gauge-equivalent to the spatial part of the previous A13 four-component raw-cochain diagnostic. Its norm RMS is checked against the previously reported \(2.3537334152693665\); the radial scalar cochain has zero ordinary face sum. The second is a genuinely different coframe and therefore defines a different finite torsion diagnostic. Its value is reported by the script/test, not identified with the previous four-dimensional result.

Run:

~~~bash
PYTHONPATH=. python scripts/audit_rgs_psu_carrier_connection.py
python -m pytest tests/test_rgs_psu_carrier_connection.py -q
~~~

## 5. Provenance and limits

The A13 generator source is preserved in the user's Library. The A13 archive manifest records the bundle SHA-256 7ca0cdfeadbe1d693d4f776a5e22d957630b206b842f4b822f3afb578c488922 and arrays R_SO3 (192,3,3), U_SU2 (192,2,2), edges (192,3). However, the original NPZ payload is not present in the repository or current Library search results. Thus this audit reproduces and compares the generator-level connection, and validates the raw-torsion RMS against the prior report; it does **not** claim bytewise verification of the archived binary matrices.

No continuum Cartan limit, TEGR identity, Lorentzian tetrad, physical time direction, Q6 tetrahedral 3-cell, or uniqueness of the PSU connection is asserted.
