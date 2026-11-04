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
unconstrained real momentum coordinates for these conventions. Factors of
two quoted with other definitions of P and q cannot be transferred unchanged.

At fixed P^2=M^2, in the centre-of-momentum frame,

    P=(M,0), q=(0,2k n), |n|=1
    k=sqrt(M^2/4-m^2).

The fixed-M two-particle mass-shell manifold has three total-momentum
coordinates on the positive mass hyperboloid and two relative-direction
coordinates on S^2. This is a five-dimensional kinematic manifold. Any
additional circle quotient or little-group restriction must be specified
separately; it does not follow from the mass-shell constraints alone.

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
