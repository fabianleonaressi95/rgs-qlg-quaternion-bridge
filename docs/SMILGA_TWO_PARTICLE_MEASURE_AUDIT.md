# Independent audit of Smilga's two-particle measure

## Scope

This implementation separates exact two-body kinematics from the additional
normalization hypothesis in Walter Smilga, *Reverse Engineering Approach to
Quantum Electrodynamics*, Journal of Modern Physics 4 (2013), 561–571,
DOI 10.4236/jmp.2013.45079.

Smilga identifies a geometrical weight with the electromagnetic coupling.
His paper explicitly describes the spin multiplicity factor 4 as heuristic
and discusses a five-dimensional integration over a basically four-dimensional
manifold. These steps remain assumptions to audit, not established results.

## Conventions and exact identities

Use eta=diag(+1,-1,-1,-1), c=hbar=1, and
P=p1+p2, q=p1-p2, so p1=(P+q)/2 and p2=(P-q)/2.
For equal mass shells p1^2=p2^2=m^2, direct algebra gives

    P.q = 0
    P^2 + q^2 = 4 m^2

The absolute Jacobian of (p1,p2) -> (P,q) is exactly 16 in eight
unconstrained real momentum coordinates. Restricting only to spatial
components gives a 6D determinant of 8. Neither raw coordinate determinant
is 2. On the fixed-invariant-mass 5D submanifold, however, an induced
Jacobian must be computed from an explicitly specified chart and induced
metric; it cannot be inferred by simply reusing either ambient determinant.
Factors of two quoted for that reduced measure therefore require their own
derivation.

At fixed P^2=M^2, in the centre-of-momentum frame,

    P=(M,0), q=(0,2k n), |n|=1
    k=sqrt(M^2/4-m^2).

The fixed-M two-particle mass-shell manifold has three total-momentum
coordinates on the positive mass hyperboloid and two relative-direction
coordinates on S^2. This is a five-dimensional kinematic manifold. Any
additional circle quotient or little-group restriction must be specified
separately; it does not follow from the mass-shell constraints alone.

## Constrained-measure reduction (new exact audit)

The full 8D Jacobian is not itself the Jacobian of the *constrained*
mass-shell measure. Define

    A = P^2 + q^2 - 4 m^2,   B = P.q
    f1 = p1^2-m^2 = (A+2B)/4
    f2 = p2^2-m^2 = (A-2B)/4.

The determinant of the map (A,B) -> (f1,f2) has absolute value 1/4.
Consequently,

    d^4p1 d^4p2 = (1/16) d^4P d^4q
    delta(f1) delta(f2) = 4 delta(A) delta(B)

and therefore

    d^4p1 d^4p2 delta(f1) delta(f2)
      = (1/4) d^4P d^4q delta(A) delta(B).

For fixed external P, the momentum-conservation delta function removes
the d^4P integral. In the CM frame P=(M,0), the remaining constraints give

    delta(P.q) = delta(M q0) = delta(q0)/M
    |q_vec| = 2k,  k=sqrt(M^2/4-m^2).

Thus the raw constrained measure is

    (1/4) d^4q delta(M^2+q^2-4m^2) delta(M q0)
      = [k/(4M)] dOmega,

and after angular integration it is pi*k/M. Multiplication by the standard
overall factor (2pi)^-2 gives

    dPhi2/dOmega = k/(16 pi^2 M),
    Phi2 = k/(4 pi M),

as expected. This gives a direct, independently testable bridge between
the 8D coordinate Jacobian, the two on-shell delta distributions and the
usual Lorentz-invariant phase space.

**Important comparison with Smilga:** the calculation above establishes
the net factor 1/4 for this explicitly stated eight-dimensional delta-shell
measure and normalization. It does not yet identify Smilga's factor 2 in
his reduced integration formula. The two factors refer to different stages
unless the intervening reduction, state normalization, angular quotient and
delta-function conventions are shown to coincide. That reconciliation
remains open; no compensating factor is inserted by hand.

With the standard Lorentz-invariant phase-space convention,

    dPhi_2=(2pi)^4 delta^4(P-p1-p2)
            product_i [d^3p_i/((2pi)^3 2E_i)]

one obtains

    dPhi_2/dOmega = k/(16 pi^2 M)
    integral_(S^2) dPhi_2 = k/(4 pi M).

These are standard kinematic results, not a derivation of alpha.

## What is computed

- Forward/inverse total-relative momentum maps.
- Exact 8D Jacobian for the declared convention.
- Equal-mass shell, orthogonality and combined-Casimir residuals.
- Centre-of-momentum configurations and nontrivial Lorentz boosts.
- Invariant two-body phase-space density and angular integral.
- Wyler's unit-radius numerical expression, kept separate.

## Open mathematical/physical steps

1. Prove that the S-matrix kernel is a genuine orthogonal projection onto the
   proposed irreducible subspace, with a specified Hilbert-space measure.
2. Prove that the bounded Lie-ball parametrization preserves the physical
   on-shell measure, including boundary and multiplicity factors.
3. Derive, rather than prescribe, the 4D-to-5D radial extension and quotient
   volume factors.
4. Derive the spin factor 4 from the actual channel sum and observable.
5. Establish the relation between the geometric weight and the renormalized
   low-energy QED coupling, including running and radiative corrections.
6. Construct a map from the finite Q6 U(1)^129 quotient to D5, if one exists.
   No such map is assumed by this audit.

## Reproducibility

Run pytest on tests/test_rgs_smilga_two_particle_audit.py and call
src.rgs_smilga_two_particle_audit.audit() for a JSON-serializable report.
Passing finite tests verifies only the exact identities listed above.


## Explicit five-coordinate chart (new)

For fixed invariant mass M>2m, introduce

    P^mu = (sqrt(M^2+|P_vec|^2), P_vec)
    n(theta,phi) = (sin(theta)cos(phi), sin(theta)sin(phi), cos(theta))
    k = sqrt(M^2/4-m^2)

and centre-of-momentum momenta

    p1* = (M/2, k n),  p2* = (M/2, -k n).

Apply the unique rotation-free boost with velocity

    v = P_vec / sqrt(M^2+|P_vec|^2)

to both momenta. This defines a five-coordinate chart
(P_x,P_y,P_z,theta,phi) -> (p1_vec,p2_vec) in R^6. It satisfies
p1^2=p2^2=m^2, P^2=M^2, P.q=0, and P_vec=p1_vec+p2_vec.
For M>2m and 0<theta<pi the chart Jacobian has rank five; the poles are
only the usual spherical-coordinate singularities. The implementation tests
this rank numerically by a central-difference Jacobian. This constructs the
kinematic five-manifold, not a D5 Lie-ball coordinate map and not a Q6 map.

The natural factorized invariant measure is

    [d^3P/(2 E_P)] dPhi_2,
    dPhi_2/dOmega = k/(16 pi^2 M),

so its density with respect to d^3P dOmega is

    k/(32 pi^3 M E_P),  E_P=sqrt(M^2+|P_vec|^2).

This is a kinematic measure identity. It does not establish that the Q6
coframe measure or Smilga's proposed projection weight equals it.
