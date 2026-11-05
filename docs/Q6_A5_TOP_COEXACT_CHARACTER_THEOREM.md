# Exact A5 representation carried by the Q6 top coexact five-mode sector

## Theorem

Let X be the six-dimensional cubical complex Q6, with cochain incidence
matrices D0 and D1 and Hodge operators

    Delta_0 = D0^T D0,
    Delta_1 = D0 D0^T + D1^T D1.

Let G=A5 be the proper icosahedral rotation group acting on the six
unoriented icosahedral axes. Its signed permutation action on the six cube
coordinates induces orthogonal cochain actions rho_0(g), rho_1(g) satisfying

    D0 rho_0(g) = rho_1(g) D0,
    D1 rho_1(g) = rho_2(g) D1.

Then the lambda=12 coexact eigenspace

    E_12^coex = ker(D0^T) intersect ker(D1^T D1 - 12 I)

has dimension five and carries exactly one copy of the irreducible
five-dimensional real representation V5 of A5.

In the standard A5 class ordering (1, 2A, 3A, 5A, 5B), its character is

    chi_E = (5, 1, -1, 0, 0),

with class sizes (1,15,20,12,12). Its character norm is

    <chi_E,chi_E> = (1/60) sum_C |C| chi_E(C)^2 = 1,

which proves irreducibility over the real/complex character theory of A5.

## Basis-independent proof by spectral projectors

The full one-cochain Hodge operator has spectrum

    Spec(Delta_1) = {2,4,6,8,10,12}.

The lambda=12 eigenspace in C^1 has dimension six: from the tensor-product
spectrum there are six choices for the unique degree-one interval factor
and the all-high degree-zero state in the other five factors.

The vertex operator has

    Spec(Delta_0) = {0,2,4,6,8,10,12},

and the lambda=12 vertex eigenspace has dimension one, since its
multiplicity is binom(6,6)=1. The map D0 is injective on this positive
eigenspace and maps it equivariantly onto the exact lambda=12 one-cochain
eigenspace. Thus

    E_12^full(C^1) = E_12^exact(C^1) orthogonal-sum E_12^coex(C^1),

and the dimensions are 6=1+5.

For a symmetric matrix A with known distinct spectrum S, its lambda
spectral projector is the Lagrange polynomial

    P_lambda(A) = product_{mu in S, mu != lambda}
                  (A - mu I)/(lambda - mu).

Because the group actions commute with the Hodge operators, these projectors
commute with the group action. The character of the coexact eigenspace is
therefore obtained without choosing any basis in the degenerate eigenspace:

    chi_coex(g)
      = tr(P_12(Delta_1) rho_1(g))
        - tr(P_12(Delta_0) rho_0(g)).

The subtraction is valid because D0 is an A5-equivariant isomorphism from
the positive-eigenvalue vertex eigenspace to the exact one-cochain
eigenspace. This formula is invariant under all choices of eigenbasis and
is an exact finite-dimensional representation identity.

The implementation constructs the two projectors as matrix polynomials,
evaluates the traces for all 60 group elements, averages within conjugacy
classes, and checks the resulting character against all five irreducible
A5 characters. The expected result is

    full lambda-12 character:   (6, 2, 0, 1, 1)
    exact lambda-12 character:  (1, 1, 1, 1, 1)
    coexact character:         (5, 1,-1, 0, 0).

The exact one-dimensional summand is trivial because A5 has no nontrivial
one-dimensional representations. Subtracting it leaves the irreducible
five-dimensional character.

## Consequence for the D5 tangent sector

The real tangent representation at the origin of the type-IV domain D5 is
the realification of its defining complex five-dimensional module. In the
RGS/Q6 construction, the five real configuration modes carry V5; their
cotangent lift is

    T*V5 = V5 direct-sum V5*.

The invariant Euclidean metric on V5 identifies V5* with V5 as A5-modules.
The induced A5 action on the ten-real-dimensional canonical phase space is
therefore the restriction of the D5 isotropy action to this embedded A5.
The complex structure pairs q and p through the oscillator phase rotation,
and A5 commutes with that SO(2) structure.

This is the precise representation-theoretic meaning of the Q6-to-D5
five-mode identification: it is not an accidental equality of dimensions
or frequencies. The finite Q6 coexact eigenspace itself is an irreducible
A5 module of the required type, and the D5 tangent sector is its canonical
cotangent/realified lift.

## Reproducibility

Run:

    pytest -q tests/test_rgs_q6_a5_character_theorem.py
    pytest -q tests/test_rgs_q6_d5_phase_space.py
    pytest -q tests/test_rgs_q6_d5_sphere_intertwiner.py

The machine-readable report is returned by
`q6_top_coexact_a5_character_report()` in
`src/rgs_q6_a5_character_theorem.py`.

## Epistemic scope

**Proved as a finite-complex representation statement:** the coexact
lambda=12 space has dimension five and character (5,1,-1,0,0), hence is
one irreducible V5 of A5.

**Numerically audited implementation:** the polynomial projectors,
60 group actions, class traces, character inner product, and decomposition
are checked with floating-point residual tolerances. The proof of the
representation statement is the projector/character argument above; the
floating-point implementation is its reproducibility audit.

**Not implied:** a global embedding of the compact Q6 gauge quotient into
the D5 bounded symmetric domain; an equality of compact Wilson and
noncompact D5 Hilbert spaces; a measure-preserving map or its Jacobian;
a physical 3+1-dimensional continuum limit; a microscopic determination of
the gauge coupling/charge normalization; or a derivation of alpha.


## Representation-theory reference

The character inner product and irreducibility criterion used here are standard finite-group representation theory; see J.-P. Serre, *Linear Representations of Finite Groups*, Graduate Texts in Mathematics 42, Springer (1977). The Q6-specific A5 action, exact spectral projector and character traces are constructed and audited in this repository. The external reference supplies the representation-theory criterion, not the Q6 matrices or their trace values.