# Oriented three-channel incidence audit

**Purpose.** This tool checks a geometrically supplied parent-to-child table against
\[
A\mapsto A+B,\qquad B\mapsto C,\qquad C\mapsto A.
\]
It does not generate or infer geometry. Passing the graph test is not evidence that
Q6/Kathara realizes the substitution.

## Input

\`oriented_children.csv\`:

\`\`\`csv
parent_id,parent_type,child_id,child_type,orientation
p001,A,c001,A,+
p001,A,c002,B,-
p002,B,c003,C,+
p003,C,c004,A,-
\`\`\`

There is one row per child occurrence. Preserve multiplicity: distinct children
of the same type must have distinct IDs and separate rows. The orientation
column is a signed incidence label; the data producer must document its geometric
convention.

Optional \`anti_twin.csv\`:

\`\`\`csv
id,antipode
p001,p101
p101,p001
\`\`\`

The antipode table must encode a fixed-point-free involution on all elements in
the incidence table. The audit checks completeness, involutivity, type
preservation, incidence invariance, and whether orientations are preserved,
reversed, or mixed.

## Exact checks

The matrix convention is rows = child type, columns = parent type, order (A,B,C):
\[
M=\begin{pmatrix}1&0&1\\1&0&0\\0&1&0\end{pmatrix}.
\]
The script checks isomorphism up to relabelling, computes the characteristic
polynomial over the integers, and enumerates graph automorphisms. The target
graph has only the identity automorphism: A is the unique type with an A-loop;
its outgoing types then distinguish B and C.

The polynomial is computed from integer trace, principal minors and determinant,
without floating point or a symbolic-algebra dependency.

## Run

\`\`\`bash
python -m pytest tests/test_oriented_three_channel_incidence.py -q
python scripts/audit_oriented_three_channel_incidence.py data/oriented_children.csv --anti-twin data/anti_twin.csv --output audit.json
\`\`\`

Exit status is 0 when the abstract substitution gate passes, 2 otherwise.
The Anti-Twin audit is reported separately; it does not change the abstract
graph-isomorphism result.

## Interpretation and limitations

- A passing report establishes that the supplied finite incidence table has the
  abstract type graph of the substitution, up to relabelling, and the
  corresponding integer matrix and characteristic polynomial.
- It does not prove that the table is derived from the RGS metric, quaternionic
  transport, Kathara geometry, or a physical refinement.
- The class-preserving Anti-Twin test is meaningful only after the antipodal
  pairing is independently defined on geometric elements.
- The script checks signed incidence labels but does not infer which geometric
  operation reverses or preserves orientation.
- The existing Q6 midpoint subdivision has transfer \(2I_3\), so it is not a
  positive instance of this test.
