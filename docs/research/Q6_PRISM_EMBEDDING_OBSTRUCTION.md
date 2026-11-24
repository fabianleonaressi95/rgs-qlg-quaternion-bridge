# Q6 inside PRISM: an embedding obstruction and the minimal product extension

**Result:** there is no injective edge-preserving graph embedding of \(Q_6\) into the existing PRISM product \(C_{64}\square C_6^3\square C_5\). This is a proved graph-theoretic obstruction, not a numerical conjecture.

## 1. The exact question

Let
\[
P=C_{64}\square C_6\square C_6\square C_6\square C_5.
\]
Does there exist an injective map \(f:V(Q_6)\to V(P)\) such that
\[
\{x,y\}\in E(Q_6)\implies \{f(x),f(y)\}\in E(P)?
\]
No requirement of inducedness is needed for the negative result. We prove that even this weaker edge-preserving embedding cannot exist.

## 2. Lemma: simple 4-cycles in products of long cycles

Every factor cycle in \(P\) has length at least five, hence has no simple 3- or 4-cycle. A simple 4-cycle in a Cartesian product of such cycles must therefore be a Cartesian rectangle: it uses exactly two distinct factor coordinates, alternating between them. In particular, the two edges incident to any one vertex of that square belong to different factor directions.

Reason: in a Cartesian product, each edge changes exactly one factor coordinate. In a closed walk of length four with four distinct vertices, every coordinate changed along the walk must return to its starting value. A single factor cannot supply a simple cycle of length four; therefore at least two coordinates are used. With only four steps and no repeated vertices, precisely two coordinates are used, each twice, in alternating order. This is the rectangle.

## 3. Theorem: no Q6 embedding into a five-factor cycle product

**Theorem.** If \(n>k\) and all cycle factors have length at least five, then \(Q_n\) has no injective edge-preserving graph embedding into a Cartesian product of \(k\) such cycles.

**Proof.** Fix a vertex \(v\in Q_n\). For every pair of distinct coordinate directions \(i,j\), the two corresponding edges incident to \(v\) belong to a 4-cycle (the \(i,j\) coordinate square). Under an injective edge-preserving map, the four vertices of this square remain distinct and its four edges form a simple 4-cycle in the target. By the lemma, the images of the two edges incident to \(f(v)\) must lie in two different Cartesian factor directions.

This applies to every pair \(i\ne j\). Thus the \(n\) edges incident to \(v\) must use pairwise distinct target factor directions. The target has only \(k\) factors, so \(n\le k\), a contradiction when \(n>k\). \(\square\)

For the existing PRISM product, \(n=6\) and \(k=5\). Hence no injective edge-preserving \(Q_6\hookrightarrow P\) exists. In particular, there is no induced Q6 subgraph and no cubical embedding of its 1-skeleton into this product.

## 4. Sharpness: six factors suffice

The factor-count obstruction is sharp. In
\[
P_6=C_{64}\square C_6^4\square C_5
\]
there are six cycle factors. Define
\[
f(v_0,\ldots,v_5)=(v_0,\ldots,v_5),
\qquad v_i\in\{0,1\}\subset \mathbb Z_{m_i},
\]
where \((m_0,\ldots,m_5)=(64,6,6,6,6,5)\).

Within each factor, the induced subgraph on \(\{0,1\}\) is a single edge. The Cartesian product of these six induced edges is exactly \(Q_6\). Thus this is an induced graph embedding. It is a result about a *modified six-factor product*, not the existing PRISM graph; adding a factor changes the model and cannot be silently treated as a coordinate relabelling.

## 5. What remains possible in the existing PRISM graph

A non-injective graph homomorphism exists. For example,
\[
g(v)=\bigl(|v|\bmod 2,0,0,0,0\bigr)\in P,
\]
where \(|v|\) is the Hamming weight and the first coordinate is in \(C_{64}\). Every Q6 edge flips parity, so its endpoints map to the adjacent first-factor vertices 0 and 1. This map has two fibers of size 32 and collapses all six cube directions onto the same edge. It is therefore a valid homomorphism but carries essentially none of Q6's six-direction geometry.

Other useful notions would need additional structure: a chain map with controlled kernel, a quotient preserving selected cycles, a correspondence of spectra, or a map from a *different* cubical complex. None follows from the parity collapse.

## 6. Reproduction

    PYTHONPATH=. python -m pytest tests/test_q6_prism_embedding_obstruction.py -q
    PYTHONPATH=. python examples/q6_prism_embedding_audit.py

The implementation exhaustively checks the explicit six-factor induced embedding and the parity homomorphism. The impossibility theorem itself is proved by the local square argument above; the code is a consistency audit, not a substitute for that proof.

## 7. Evidence ledger

| Statement | Status |
|---|---|
| Every simple 4-cycle in a Cartesian product of cycles of length \(\ge5\) is a rectangle in two factors | PROVED |
| \(Q_n\) cannot embed injectively and edge-preservingly into a product of \(k<n\) such cycles | PROVED |
| Existing PRISM \(C_{64}\square C_6^3\square C_5\) has no Q6 embedding | PROVED |
| Six cycle factors suffice; the binary corner subgraph is induced Q6 | PROVED |
| Parity collapse is a graph homomorphism with two 32-vertex fibers | PROVED |
| A useful Q6-to-existing-PRISM chain/spectral/transport map exists | OPEN |

No physical, Green–Kubo, or dynamical conclusion is asserted.
