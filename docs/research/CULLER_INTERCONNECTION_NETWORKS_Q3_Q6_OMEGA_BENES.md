# Interconnection networks: Culler, Q3/Q6, Omega, Beneš and Clos

**Research recovery record**  
**Recovered session:** 2026-11-16  
**Repository:** `rgs-qlg-quaternion-bridge`  
**Status:** bibliography recovered; historical finite-network results recorded; independent rerun of the Omega/Beneš enumeration remains pending.

## 1. Starting point: David E. Culler

The primary lecture source is:

- David E. Culler, **“Interconnection Network Topology Design Trade-offs”**, UC Berkeley, CS 258, Spring 1999.
- [Original Berkeley PDF — lecture handout](https://people.eecs.berkeley.edu/~culler/cs258-s99/lec18.pdf)
- [Slide transcription / presentation record](https://slidetodoc.com/interconnection-network-topology-design-tradeoffs-cs-258-spring/)

The lecture compares linear arrays and rings, multidimensional meshes and tori, trees and fat-trees, hypercubes and crossbars. It uses graph and implementation metrics including degree, diameter, average distance, bisection bandwidth, wire length, routing delay, latency and throughput. It also discusses dimension-order routing, contention, channel width and physical pin-count constraints.

A complementary lecture handout is:

- David E. Culler, **“Scalable Synchronization Operations / Scalable, High Performance Interconnection Networks”**, UC Berkeley CS 258, Spring 1999.
- [Berkeley PDF — lecture handout](https://people.eecs.berkeley.edu/~culler/cs258-s99/lec30.pdf)

The latter includes a topology comparison table and discussion of routing restrictions, deadlock freedom, channel-dependence graphs, deterministic versus adaptive routing, switch design and flow control.

## 2. Culler–Singh–Gupta textbook

- David E. Culler, Jaswinder Pal Singh and Anoop Gupta, *Parallel Computer Architecture: A Hardware/Software Approach*, Morgan Kaufmann, 1998.
- [Chapter 10 — Interconnection Network Design](https://www.oreilly.com/library/view/parallel-computer-architecture/9780080573076/chapter-192.html)

Chapter 10 covers basic definitions, communication performance, organizational structure, interconnection topologies, design trade-offs and routing. This is the natural extended reference for the lecture material. The chapter page is a publisher/platform chapter record, not a freely downloadable full-book PDF.

## 3. Foundational multistage-network references

### Clos networks

- Charles Clos, “A Study of Non-Blocking Switching Networks,” *Bell System Technical Journal* 32(2), 406–424 (1953).
- DOI: [10.1002/j.1538-7305.1953.tb01433.x](https://doi.org/10.1002/j.1538-7305.1953.tb01433.x)
- [Bell Labs/Nokia publication record](https://www.nokia.com/bell-labs/publications-and-media/publications/a-study-of-non-blocking-switching-networks/)
- [Scanned complete BSTJ volume, March 1953 (article begins p. 406)](https://bitsavers.trailing-edge.com/magazines/Bell_System_Technical_Journal/BSTJ_V32N02_195303.pdf)

Clos's construction is the classical reference for strictly nonblocking multistage switching under its stated switch-size and stage conditions.

### Beneš networks

- Václav E. Beneš, “Algebraic and Topological Properties of Connecting Networks,” *Bell System Technical Journal* 41(4), 1249–1274 (1962).
- DOI: [10.1002/j.1538-7305.1962.tb03277.x](https://doi.org/10.1002/j.1538-7305.1962.tb03277.x)
- Václav E. Beneš, “Optimal Rearrangeable Multistage Connecting Networks,” *Bell System Technical Journal* 43(4), 1641–1656 (1964).
- DOI: [10.1002/j.1538-7305.1964.tb04103.x](https://doi.org/10.1002/j.1538-7305.1964.tb04103.x)

The 1964 paper is the primary reference for rearrangeable multistage networks. Rearrangeability means that any input-output permutation can be realized after rearranging existing paths; it is not the same property as strict nonblocking.

### Survey and routing context

- “Multistage Interconnection Networks for Multiprocessor Systems,” *Advances in Computers* 26 (1987), 155–199.
- DOI: [10.1016/S0065-2458(08)60007-8](https://doi.org/10.1016/S0065-2458(08)60007-8)
- [Publisher record and abstract](https://www.sciencedirect.com/science/article/pii/S0065245808600078)

The survey distinguishes topology, switch functionality and control strategy, and discusses blocking, nonblocking, rearrangeable and multipath networks, including Omega routing and conflicts.

- Kuo-Yu Chen, *Multistage interconnection networks: improved routing algorithms and fault tolerance*, M.S. thesis, New Jersey Institute of Technology (1991).
- [Institutional thesis record](https://digitalcommons.njit.edu/theses/1291/)

The thesis record describes routing and fault-tolerance work on Beneš and Clos networks, including permutation representations and edge-colouring methods.

### Modern engineering reference

- William J. Dally and Brian Towles, *Principles and Practices of Interconnection Networks*, Morgan Kaufmann (2004).
- [Elsevier book record](https://www.educate.elsevier.com/book/details/9780122007514)
- [Contents and chapter access](https://www.oreilly.com/library/view/principles-and-practices/9780122007514/)

Relevant chapters include topology basics, butterfly and torus networks, nonblocking networks, routing, adaptive routing, flow control, deadlock/livelock and performance analysis. This is a book reference, not an open PDF.

## 4. Mathematical bridge to Q3 and Q6

Keep two mathematical objects separate:

1. **Hypercube graph:** a graph whose vertices are bit strings and whose edges join strings at Hamming distance one.
2. **Multistage switching network:** a directed, layered arrangement of switching elements that realizes a specified family of input-output permutations.

They are related by binary addressing and routing ideas, but they are not the same graph.

### Q3

The three-cube graph has
[
V(Q_3)=\mathbb F_2^3,\qquad |V|=8,\qquad |E|=12,
\]
and degree 3. Its vertices can be indexed by the eight I Ching trigram bit strings.

### Q6

The six-cube graph has
[
V(Q_6)=\mathbb F_2^6,\qquad |V|=64,\qquad |E|=192,
\]
and degree 6. With the product decomposition
[
Q_6\cong Q_3\square Q_3,
\]
a six-bit vertex is an ordered pair of three-bit strings: lower and upper trigram under the repository's declared bit-order convention.

Its graph distance is Hamming distance:
[
d(x,y)=\sum_{j=0}^{5}|x_j-y_j|.
\]
A dimension-order route flips the differing coordinates in a fixed coordinate order. Different coordinate orders may yield different shortest paths, so a routing convention must be declared.

The exact graph-Laplacian spectrum is
[
\lambda_k=2k,\qquad
m_k=\binom 6k,\qquad k=0,\ldots,6.
\]
This is a finite graph identity, not a statement about the performance of a physical interconnection network.

## 5. Recovered Omega and Beneš 8×8 experiment

The 2026-11-16 research session recorded the following architecture counts:

| Network | Stages | 2×2 switches per stage | Total 2×2 switches | Routing property |
|---|---:|---:|---:|---|
| Omega 8×8 | 3 | 4 | 12 | Blocking |
| Beneš 8×8 | 5 | 4 | 20 | Rearrangeable |

The session also recorded exhaustive permutation-enumeration outputs:

| Check | Recorded result | Interpretation |
|---|---:|---|
| All 8×8 permutations | 8! = 40,320 | Enumeration domain |
| Omega routable permutations | 4,096 / 40,320 | Output of the historical routing implementation |
| Beneš routable permutations | 40,320 / 40,320 | Output of the historical routing implementation |
| Bit-reversal permutation on Omega | Not routable | Historical test result |
| Bit-reversal permutation on Beneš | Routable | Historical test result |

**Evidence status:** these are recovered session results, not newly executed tests in this commit. The exact historical Python files and raw logs have not yet been located in this repository. The enumeration must be independently reimplemented or recovered, then rerun with explicit switch numbering, interstage permutation, route-setting algorithm and conflict semantics before these counts are promoted to reproducible repository evidence.

In particular, do not infer that every 8×8 Omega convention has exactly the same routable-permutation count without fixing its wiring convention. The Omega count is topology- and convention-dependent.

## 6. Connection to the existing Q6–Green–Kubo–supergolden dossier

The finite Q6 / Green–Kubo / supergolden implementation is already present:

- [Mathematical evidence record](../Q6_HYPERCUBE_GREEN_KUBO_SUPERGOLDEN.md)
- [Implementation](../../src/q6_hypercube_gk_supergolden.py)
- [Tests](../../tests/test_q6_hypercube_gk_supergolden.py)
- [Reproducible audit example](../../examples/q6_hypercube_gk_supergolden_audit.py)

That dossier contains finite Q6 cochain and cosine identities, a closed-2-cochain projector, finite-dimensional OU Green–Kubo formulas and exact algebra concerning the supergolden root. It explicitly leaves physical PSU viscosity, continuum Yang–Mills convergence and universal supergolden scaling open.

The network-topology literature supplies the architecture/routing side; the Q6 dossier supplies the finite cubical geometry and conditional transport side. A mathematical bridge between them still requires a defined map between the cubical graph/cochain model and the switching-network state/routing model. No such equivalence is asserted here.

## 7. Reproduction and next verification gate

Existing Q6 dossier tests:
```bash
PYTHONPATH=. python -m pytest tests/test_q6_hypercube_gk_supergolden.py -q
PYTHONPATH=. python examples/q6_hypercube_gk_supergolden_audit.py
```

Required new independent network tests:
1. Generate the exact Omega 8×8 and Beneš 8×8 stage wiring from a written convention.
2. Enumerate all 40,320 permutations.
3. Implement a route-setting algorithm independently of the historical implementation.
4. Verify no switch is assigned conflicting outputs and every input/output is used exactly once.
5. Check bit reversal and additional adversarial permutations.
6. Record raw JSON output, source revision and algorithm details.
7. Compare the Q3/Q6 Hamming-routing model with the multistage network only through an explicitly defined routing map.

## 8. Evidence policy

- **PROVED:** finite combinatorial identities with a mathematical derivation.
- **NUMERICALLY VERIFIED:** results reproduced by committed code and tests.
- **RECOVERED / NOT RERUN:** historical outputs whose original execution record is not yet available.
- **OPEN:** proposed cross-model equivalence, continuum limit or physical interpretation.

No copyrighted book or publisher PDF is copied into this repository. The bibliography records stable DOI, publisher, institutional and author-hosted links; access restrictions and source type are stated explicitly.
