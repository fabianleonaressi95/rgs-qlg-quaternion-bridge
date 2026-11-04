# Q6 Wilson microhistories, symmetric quotient, and local modal semantics

## Exact combinatorics

Fix an oriented edge set E with |E|=192 and an integer flux
n=(n_e) in Z^E. Put s=||n||_1=sum_e |n_e| and m_e=|n_e|.
A labelled microhistory is a word of length s whose signed token
sgn(n_e)e occurs exactly m_e times. The symmetric group S_s acts by
permuting the s slots. Its stabilizer is the product of the symmetric
groups S_(m_e), hence has order product_e m_e!. Orbit-stabilizer gives

    |S_s . h| = s! / product_e m_e!,
    |S_s . h| / s! = 1 / product_e m_e!.

Consequently the factorial Bessel coefficient majorant from PR #156 is
exactly expressible as

    (beta/2)^s / product_e |n_e|!
      = (beta/2)^s * |S_s . h| / s!.

This is an exact orbit-stabilizer identity, not a heuristic analogy.

## What this does and does not say about the Wilson measure

The identity concerns the factorial majorant of the Bessel product. It does
not show that the exact product of modified Bessel coefficients is itself
the pushforward of uniform counting measure on labelled histories. The
Bessel series also contains pair-creation indices j_e with weights
(beta/2)^(2 j_e)/(j_e! (m_e+j_e)!). A measure-level identification must
include these indices and their weights. The orbit formula identifies the
factorial denominator in the majorant; it does not identify the full Wilson
measure with a quotient counting measure.

Different signed occupancy vectors are distinct quotient orbits. Permuting
event labels removes only the arbitrary ordering of the unit-flux events.

## Local Kripke semantics

For a fixed history h, take the worlds to be its S_s orbit and define

    w R_h v iff v = sigma.w for some sigma in S_s.

This is an equivalence relation. The induced Kripke frame is an S5 frame.
For any predicate P,

    w |= Box_(S_s) P iff for every v in S_s.h, P(v).

Thus, on this frame, necessity means invariance throughout the relabelling
orbit. The code checks the equivalence-frame property and invariant versus
non-invariant predicates exactly.

## Separation from Gödel's ontological proof

This local Box_(S_s) is not automatically Gödel's Box. Gödel's ontological
argument uses a modal language with positivity, essence, and necessary
existence predicates, together with specific modal axioms (commonly
presented in an S5 setting). A permutation orbit supplies an S5 accessibility
relation, but supplies none of those positivity axioms or an existence
predicate. Orbit-stabilizer does not prove that any object exists necessarily,
nor does it prove consistency of a positive-property system.

The defensible bridge is conditional and structural:

1. **Proved:** event relabellings form S_s; orbit-stabilizer gives
   s!/product_e m_e!.
2. **Proved:** the factorial majorant is
   (beta/2)^s |S_s.h|/s!.
3. **Proved:** orbit accessibility is an equivalence relation and yields
   local S5 semantics.
4. **Open:** exact Wilson Bessel measure equals a quotient counting measure;
   Bessel pair-creation labels and weights must be incorporated.
5. **Not proved:** Gödel's positivity/essence/necessary-existence axioms
   follow from Q6, RGS, or Wilson dynamics.

## Reproducibility

Run: python scripts/run_q6_wilson_microhistory_modal_bridge.py

The JSON report uses a small four-edge exact example. It is a structural
unit test, not a truncation or numerical evaluation of the 124-dimensional
Q6-D5 fiber.
