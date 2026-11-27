# GEC-6 and the icosahedral matching orbits

## Scope and conventions

This dossier tests whether the lower/upper trigram split and its diagonal
matching are selected by the six unoriented icosahedral axes already used by
the signed-axis intertwiner. Coordinates are those of
`build_icosahedral_c6()`; Q6 bit ordering remains the repository's declared
big-endian convention.

There are three different actions and they must not be conflated:

1. the permutation projection of the 60 signed-axis matrices on the six
   *unoriented axis labels*;
2. the full signed action on Q6 vertices, where a negative sign complements
   the corresponding bit;
3. the induced action on the GEC-6 quotient
   `c: Λ²(Z⁶) → A5`.

The orbit counts below are exact finite enumeration once the 60 signed
permutation matrices have been reconstructed. Their extraction from floating
icosahedral coordinates is numerically certified by the existing intertwiner
audit; this dossier does not reclassify that floating reconstruction as exact
arithmetic.

## A5 orbit structure on the six axis labels

Let (p_gin S_6) be the permutation part of the signed-axis action. Its
action on unordered coordinate pairs, 3+3 partitions and perfect matchings
gives:

| Object | A5 orbit sizes |
|---|---:|
| Unordered pairs of axes (15) | 15 |
| Unordered partitions into two triples (10) | 10 |
| Perfect matchings (15) | 5 + 10 |

Thus A5 is transitive on all 15 pair directions and on all ten 3+3
partitions. No individual lower/upper partition is A5-invariant.

The 15 perfect matchings split into one distinguished orbit of five and a
generic orbit of ten. The stabilizer of a matching in the five-element orbit
has order (60/5=12) (an (A_4) subgroup); a generic matching has stabilizer
of order (60/10=6).

## The existing C6 tensor detects the distinguished five-orbit

In the deterministic axis orientation returned by `build_icosahedral_c6()`,
the C6 tensor has integral off-diagonal entries (C_{ij}in{+1,-1}).
For a perfect matching (M={{i_1,j_1},{i_2,j_2},{i_3,j_3}}), define
[
chi_{C_6}(M)=C_{i_1j_1}C_{i_2j_2}C_{i_3j_3}.
]

The finite audit gives:

- (chi_{C_6}(M)=+1) on exactly the five-element A5 orbit;
- (chi_{C_6}(M)=-1) on the ten-element orbit.

Changing the sign convention of one or more oriented axis representatives
multiplies every matching product by the same global factor
(prod_ieta_i). It may exchange the labels (+1) and (-1), but it does
not change the 5+10 orbit partition. The orbit distinction is intrinsic to
the A5 permutation action; the sign assigned to each class depends on the
chosen axis representatives.

The diagonal matching
[
M_0={{0,3},{1,4},{2,5}}
]
is in the distinguished five-orbit. Hence it is **icosahedrally special**, not
a generic matching. But it is not uniquely selected: for every 3+3 partition
there are exactly **two** special perfect matchings whose three pairs cross
between the two blocks. For the repository's lower/upper partition these are
[
{{0,3},{1,4},{2,5}},qquad
{{0,4},{1,5},{2,3}}.
]
The full A5 stabilizer of the lower/upper partition has order six and
interchanges these two special cross-matchings. Therefore the icosahedral
data plus the unordered lower/upper split do not choose one of the two.

This is the positive geometric result: the existing C6 tensor identifies a
canonical *family of five* matchings, and the chosen diagonal is a member of
that family. It is not a unique lower-to-upper identification.

## Compatibility with the GEC-6 quotient: an important obstruction

The signed-axis action on Q6 is not automatically an action on the GEC-6
quotient. Write a signed permutation as
[
g e_i=epsilon_i e_{p(i)},qquad epsilon_iin{pm1}.
]
For (i<j<k), the triangle relation
[
r_{ijk}=e_iwedge e_j+e_jwedge e_k-e_iwedge e_k
]
lies in (ker c), where (c(e_iwedge e_j)=e_j-e_i). Applying the signed
permutation and contracting gives coefficients proportional to
[
epsilon_i(epsilon_k-epsilon_j),quad
epsilon_j(epsilon_i-epsilon_k),quad
epsilon_k(epsilon_j-epsilon_i)
]
on the three distinct target coordinates. This vanishes for every triple
if and only if all six signs (epsilon_i) are equal. Hence the signed action
preserves (ker c) exactly when its signs are all (+1) or all (-1).

Among the 60 reconstructed icosahedral signed-axis actions, exactly six
satisfy this criterion. They are exactly the stabilizer of the selected
3+3 lower/upper partition. The quotient-compatible subgroup still moves
(M_0) through two special cross-matchings; it does not select one uniquely.

Consequently:

- the full A5 orbit calculation on matchings is a statement about the
  permutation action on the six *unoriented axis labels*;
- the full signed Q6 action does **not** descend to the GEC-6 contraction
  quotient;
- the six-element quotient-compatible subgroup preserves the block
  partition but not the individual diagonal matching.

This distinction prevents an invalid claim of full A5-equivariance for the
GEC-6 quotient. The existing GEC-6 contraction and the signed Q6
icosahedral action remain different structures until an alternative
A5-equivariant quotient map is constructed.

## Evidence status

| Claim | Status |
|---|---|
| A5 pair, partition and matching orbit sizes (15; 10; 5+10) | FINITE-EXACT after group extraction |
| C6 matching-product invariant separates the 5- and 10-orbits | FINITE-EXACT in the declared axis convention |
| Diagonal matching belongs to the distinguished five-orbit | FINITE-EXACT |
| Exactly two special cross-matchings per 3+3 partition | FINITE-EXACT |
| Signed action preserves (ker c) iff all six signs are equal | PROVED algebraically |
| Exactly six of the 60 signed actions descend to the GEC-6 quotient | FINITE-EXACT |
| A unique lower-to-upper matching selected by A5 and the block split | DISPROVED |
| Physical interpretation of the two special cross-matchings | OPEN |

## Reproduction

```bash
PYTHONPATH=. python -m pytest tests/test_gec6_icosahedral_matching_orbits.py -q
PYTHONPATH=. python -m src.gec6_icosahedral_matching_orbits
```
