# Unified PSU–RGS–Anti-Twin wormhole: implementation crosswalk

**Date:** 2026-11-28  
**Primary mathematical specification:** [theorem repository, Unified PSU–RGS–Anti-Twin wormhole construction](https://github.com/fabianleonaressi95/theorem/blob/rgs-q6-antitwin-spacetime-bridge/research/UNIFIED_PSU_RGS_ANTITWIN_WORMHOLE_CONSTRUCTION_2026-11-28.md).

This companion note makes the existing finite Cartan implementation serve the single PSU two-sheet/wormhole programme. It does not claim that a Q6-derived Lorentzian wormhole or its field equations have been implemented.

## 1. Shared data model

The intended finite object is an oriented four-dimensional refinement \(K_h\) with 4-cells and complete attaching maps; two sheet labels \(s\in\{+,-\}\), exchanged by \(J_{AT}=\tau_x\); a throat interface subcomplex \(\Sigma_h\); independent Lorentz links \(L_{vw}\in SO^+(1,3)\); oriented edge coframes \(E_{vw}\in\mathbb R^{1,3}\); and spin lifts of the same links only after an oriented/time-oriented spin structure is available.

The current Q6 cubical data provide exact six-cube counts \((64,192,240,160,60,12,1)\). The 4-skeleton is not a 4-manifold by itself. A future constructor must output selected 4-cells, all incidence maps, local links, throat gluing and a refinement family.

## 2. Reuse existing Cartan core

Reuse, without changing their established scope:
- scripts/discrete_cartan_core.py: Lorentz-link and affine edge operations, face holonomy, transported coframe closure, gauge transformations, affine products/inverses and tetrahedral boundary audit;
- tests/test_discrete_cartan_core.py: finite gauge covariance, independence of face curvature and closure defect, exact affine boundary identity and 4-simplex facet checks;
- scripts/audit_cartan_refinement_consistency.py and tests/test_cartan_refinement_consistency.py: manufactured smooth-sampling refinement audit.

The present code checks finite affine Cartan algebra. It does not yet implement a Q6 4D manifold, two-sheet throat gluing, a variational action, physical energy conditions or a spinor operator.

## 3. Required extensions, in dependency order

### WQ6-01 — Q6 four-dimensional resolution
Input: full cubical boundary operators of the six-cube. Output: a 4D refinement with verified local links, orientations and shape regularity. Reject if a codimension-one face has the wrong incidence or if an interior vertex link is not the required 3-sphere (boundary links are 3-balls), unless a separately defined controlled singular stratum is intended.

### WQ6-02 — two-sheet throat gluing
Input: WQ6-01 and the PSU/Anti-Twin sheet involution. Output: throat subcomplex, transition maps and invertible gluing. Check inverse gluing, orientation/time-orientation compatibility, and consistency of coframe and connection pullbacks.

### WQ6-03 — coframe regularity
Input: throat-adjacent 4-cells and their Lorentz frames. Output: coframe on each 4-cell and induced metric. Check rank four, signature \((+---)\), a lower bound on the smallest singular value in a fixed normalized frame, finite metric/curvature/torsion observables and chart compatibility.

### WQ6-04 — independent Cartan connection
Input: coframe and gluing. Output: independent \(SO^+(1,3)\) links and affine edge data. Reuse existing face_holonomy and face_torsion conventions; retain the exact affine tetrahedral boundary identity. Test that curvature and torsion can vary independently.

### WQ6-05 — one action, two independent variations
Implement the discrete Palatini–Cartan face action with explicit branch-controlled matrix logarithms, transported area bivectors, cosmological term, matter term and throat/interface term. Derive and implement separate gradients with respect to coframe and Lorentz links. Validate analytic gradients against finite differences and gauge covariance. Do not replace the derivative of log(H) by H^{-1}dH for noncommuting matrices.

### WQ6-06 — wormhole and limiting solutions
Implement the earlier toroidal throat ansatz as the primary historical branch:
\[
ds^2=-N(l)^2c^2dt^2+dl^2+r(l)^2d\theta^2+
[R(l)+r(l)\cos\theta]^2d\varphi^2.
\]
Keep the spherical two-ended ansatz as a separate comparison. Check regularity, no horizon on a traversable branch, tidal forces, energy conditions from the same action, and throat flux balance. Then construct a controlled decoupling family and compare its limit with Schwarzschild in regular Kruskal coordinates. The Einstein–Rosen limit is non-traversable and must not be labelled an open traversable wormhole.

### WQ6-07 — spinorial throat transport
Lift the same Lorentz links to \(Spin^+(1,3)\), define two-sheet spinor matching, and construct a refinement-dependent Dirac operator. Audit locality, spectral/resolvent convergence, doubling, chirality, anomalies and regularity across the throat. T149 rules out simply reusing the bare Q6 Clifford operator as the desired massless Dirac operator.

### WQ6-08 — convergence of solutions
Only after WQ6-01–07: prove uniform stability/compactness modulo gauge and identify limits of discrete critical points. Existing smooth-sampling consistency is a local consistency theorem, not this solution-convergence result.

## 4. Toroidal throat data and cautions

The earlier PSU–Anti-Twin model used a toroidal throat, with PSU and Anti-Twin on opposite sides, a proposed one-way current, and external Vesica-Piscis anisotropy triplets. The one-way current and its orientation remain hypotheses until derived from the discrete/continuum equations. The sheet involution alone cannot produce a preferred direction.

For the toroidal chart require \(N>0\), \(r>0\), and \(R-r>0\), with smooth matching and finite invariants. The spherical comparison \(R(u)=\sqrt{u^2+r_0^2}\) is not a replacement for the toroidal geometry.

A concrete local coframe is
\[
e^0=N(l)c\,dt,\quad e^1=dl,\quad
e^2=r(l)d\theta,\quad e^3=A(l,\theta)d\varphi,\quad
A=R(l)+r(l)\cos\theta.
\]
Its volume density is \(NcrA\). Thus \(N,r,A>0\) gives rank four, and \(R-r>0\) is a sufficient global angular-patch condition. For even throat profiles \(N,r,R\), the coframe remains nondegenerate at \(l=0\). A reference torsion-free connection has
\[
\omega_{\rm LC}{}^0{}_1=(N'/N)e^0,\quad
\omega_{\rm LC}{}^2{}_1=(r'/r)e^2,\quad
\omega_{\rm LC}{}^3{}_1=(\partial_lA/A)e^3,\quad
\omega_{\rm LC}{}^3{}_2=(\partial_\theta A/(rA))e^3.
\]
The independent connection is \(\omega=\omega_{\rm LC}+K\), \(K_{ab}=-K_{ba}\), so
\[
T^a=K^a{}_b\wedge e^b,\qquad
R(\omega)=R(\omega_{\rm LC})+D_{\rm LC}K+K\wedge K.
\]
This reference decomposition does not impose \(K=0\); independent variation remains mandatory.

A minimal throat-localized test connection is \(K^{02}=\kappa(l)\,dl\), \(K^{20}=-K^{02}\). It gives nonzero torsion components \(T^0=-\kappa\,dl\wedge e^2\) and \(T^2=-\kappa\,dl\wedge e^0\) (for \(\eta=(+---)\) and the stated internal-index convention). Along a radial path the reference connection has no \(dl\) component, so the test transport is
\[
\mathcal P_\gamma=\exp\!\left[-J_{02}\int_\gamma\kappa(l)\,dl\right].
\]
This is a kinematic manufactured ansatz for testing throat transport and finite torsion, not a solution of the action. Its curvature, interface matching and field equation for \(\kappa\) still have to be derived.

## 5. Status boundary

**Already available:** exact Q6 cubical counts; T129/T161 golden frame; T148 two-sheet algebra; T149 bare-operator obstruction; finite affine Cartan core and tests; exact affine tetrahedral boundary identity; local smooth-sampling Cartan consistency proof and manufactured RK4 audit.

**Still open:** a Q6-derived manifold-like 4D refinement; Lorentzian coframe/time selector; throat gluing; full discrete action and independent variations; an actual traversable solution and its energy analysis; Schwarzschild–Kruskal degeneration theorem; Spin+(1,3) spinor construction; convergence of discrete solutions.

The implementation and theorem graph must point to this same construction. A family/topic edge is not automatically a proof dependency.
