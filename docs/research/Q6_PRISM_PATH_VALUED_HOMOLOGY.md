# A path-valued Q6 graph map surjecting onto PRISM periodic H1

**Status: PROVED existence construction.** The construction deliberately relaxes the edge-to-edge condition: a source edge may map to a finite target edge path. It uses Q6 as a graph only, not as the filled cubical 6-cube.

## 1. Source and target

Let \(G=Q_6^{(1)}\), with 64 vertices and 192 edges. Its integral graph cycle group is free abelian of rank
\[
\beta_1(G)=192-64+1=129.
\]
Let \(P=C_{64}\square C_6^3\square C_5\) with its standard cubical structure. Its periodic first homology is
\[
H_1(P;\mathbb Z)\cong\mathbb Z^5.
\]
Choose the common basepoint \(p=(0,0,0,0,0)\). Each factor gives a coordinate loop \(\gamma_i\) of length \(n_i\in(64,6,6,6,5)\), representing the \(i\)-th standard generator of \(\mathbb Z^5\).

## 2. A deterministic fundamental cycle basis of Q6

Identify Q6 vertices with integers 0 through 63. For \(v\ne0\), define its parent by deleting its highest nonzero bit:
\[
\operatorname{par}(v)=v\mathbin{\mathrm{xor}}2^{\lfloor\log_2v\rfloor}.
\]
These 63 edges form a spanning tree \(T\). The remaining 129 edges are chords \(c_0,\ldots,c_{128}\), ordered lexicographically by their endpoint pairs.

The fundamental cycles associated with these chords form an integral basis of the graph cycle group. For each chord \(c_j=(u,v)\), its fundamental cycle consists of that chord and the unique tree path returning from \(v\) to \(u\), with orientation fixed by the ordered chord.

## 3. Define the path-valued map

Send every source vertex to \(p\). For the first five chords, assign
\[
c_j\longmapsto\gamma_j,\qquad j=0,\ldots,4.
\]
Send every tree edge and every remaining chord to the constant path at \(p\). Since every edge path begins and ends at \(p\), these assignments extend continuously over the entire source graph.

The induced map on the fundamental cycle basis is exactly
\[
c_j^*\longmapsto
\begin{cases}
e_j\in\mathbb Z^5,&0\le j<5,\\
0,&5\le j<129.
\end{cases}
\]
Consequently the induced homomorphism
\[
\boxed{H_1(G;\mathbb Z)\longrightarrow H_1(P;\mathbb Z)=\mathbb Z^5}
\]
is represented by the \(5\times129\) matrix \([I_5\;0]\), and is surjective of rank five.

This is an explicit realization of nonzero periodic winding after replacing edge-to-edge maps by edge paths.

## 4. What this construction proves — and does not prove

It proves that the four-step obstruction is exactly tied to the single-edge image condition: allowing edge paths permits a surjective map on periodic \(H_1\). The five selected source edges are replaced by paths of lengths 64, 6, 6, 6 and 5, respectively; all other edges are collapsed to the basepoint.

It is intentionally a universal, basis-dependent existence construction. It does **not** preserve Q6 distances, adjacency locally, the I-Ching line/trigram interpretation, the cube's symmetries, or any Q6 dynamics. The choice of five chords is deterministic but not canonical from the I-Ching structure. It must not be interpreted as a physical embedding or a transport law.

## 5. Why the 240 plaquettes cannot all be filled

This is a map from the Q6 graph \(G\), not from the filled cubical 6-cube \(Q_6\). If it extended continuously over all 240 square faces, its induced map on \(H_1\) would factor through
\[
H_1(Q_6;\mathbb Z)=0,
\]
and hence would be zero. Since the constructed graph map is surjective onto \(\mathbb Z^5\), such an extension is impossible. At least some source plaquette boundaries map to non-null-homologous loops and cannot be filled in the PRISM target.

This identifies the modelling decision required for further work: which Q6 plaquettes are retained as filled 2-cells, and which cycles are intended to remain as graph-level generators?

## 6. Reproduction

    PYTHONPATH=. python -m pytest tests/test_q6_prism_path_valued_homology.py -q
    PYTHONPATH=. python examples/q6_prism_path_valued_homology_audit.py

## 7. Evidence ledger

| Claim | Status |
|---|---|
| Q6 graph cycle rank is 129 | PROVED |
| The parent rule defines a spanning tree | PROVED |
| The 129 non-tree edges give a fundamental cycle basis | PROVED |
| The path-valued map is continuous on the graph | PROVED by construction |
| Induced \(H_1\) matrix is \([I_5\;0]\) | PROVED |
| Induced map is surjective onto \(\mathbb Z^5\) | PROVED |
| Extension over all 240 filled plaquettes | IMPOSSIBLE |
| Canonical, symmetry-compatible, or dynamical interpretation | OPEN |
