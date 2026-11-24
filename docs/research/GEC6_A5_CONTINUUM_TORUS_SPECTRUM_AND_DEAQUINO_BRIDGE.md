# GEC-6: continuum A5 torus spectrum and De Aquino threshold

**Status:** the flat-torus spectrum and its first eigenspace are proved. The comparison with De Aquino is a conjectural normalization, not a physical identification.

## 1. Flat torus and dual spectrum

Let
\[
V=\{x\in\mathbb R^6:\sum_{i=1}^6x_i=0\},\quad
A_5=V\cap\mathbb Z^6,\quad \mathbb T_{A_5}=V/A_5.
\]
The Euclidean metric descends to a compact flat five-torus. Its character lattice is the metric dual lattice A5*. For \(\xi\in A_5^*\), the nonnegative Laplacian has eigenvalue
\[
\lambda(\xi)=4\pi^2\|\xi\|^2.
\]

In the simple-root basis, let C be the A5 Cartan matrix. The dual Gram matrix is C^{-1}; setting Q=6C^{-1} gives
\[
Q=\begin{pmatrix}
5&4&3&2&1\\4&8&6&4&2\\3&6&9&6&3\\
2&4&6&8&4\\1&2&3&4&5
\end{pmatrix}.
\]
Hence
\[
\|\xi_k\|^2=\frac16 k^TQk,\qquad
\boxed{\lambda_k=\frac{2\pi^2}{3}k^TQk}.
\]

## 2. Exact first eigenvalue and multiplicity

The minimum of k^TQk over nonzero integer k is 5. The leading principal minors of Q-(3/2)I are
\[
7/2,\quad27/4,\quad81/8,\quad189/16,\quad243/32,
\]
so Q is strictly greater than (3/2)I. If k^TQk<=5, then ||k||^2<10/3, so every integer coordinate lies in {-1,0,1}. The exhaustive finite set {-1,0,1}^5 has 243 elements; exact integer evaluation gives minimum 5 and exactly twelve minimizers:
\[
\begin{aligned}
&(-1,0,0,0,0),\;(-1,1,0,0,0),\;(0,-1,1,0,0),\\
&(0,0,-1,1,0),\;(0,0,0,-1,1),\;(0,0,0,0,-1),
\end{aligned}
\]
and their six negatives. Therefore
\[
\boxed{\lambda_1(\mathbb T_{A_5})=10\pi^2/3,
\qquad \operatorname{mult}(\lambda_1)=12.}
\]
This multiplicity is an exact shortest-vector count, not a fitted numerical degeneracy.

## 3. Discrete-to-continuum relation

For the normalized Laplacian of the finite A5-root Cayley quotients, the general-modulus result is
\[
\lim_{m\to\infty}m^2\gamma_m=2\pi^2/3.
\]
With the continuum normalization above,
\[
\boxed{\lim_{m\to\infty}m^2\gamma_m
=\frac15\lambda_1(\mathbb T_{A_5}).}
\]
The factor 1/5 comes from averaging the 30 root directions as 15 opposite pairs. This is a spectral scaling relation, not a physical field-theory claim.

## 4. De Aquino: exact threshold and proposed bridge

The mass-ratio expression used in the earlier project audit is
\[
\frac{m_g}{m_{i0}}=3-2\sqrt{1+\xi^2}.
\]
Its zero is exactly
\[
3=2\sqrt{1+\xi^2}
\quad\Longleftrightarrow\quad
\boxed{\xi_*^2=5/4,\;|\xi_*|=\sqrt5/2}.
\]
The independent A5 result is min_{k nonzero} k^TQk=5. The two values can be formally aligned by defining
\[
\xi_k^2:=\frac14 k^TQk.
\]
For the twelve shortest dual vectors this gives xi_k^2=5/4, the De Aquino zero threshold.

**This is a candidate bridge only.** The A5 torus does not determine the physical momentum scale, the units in xi=Delta p/(m_i0 c), or the factor 1/4. No canonical map from physical momentum to A5*, no operator intertwining the torus Laplacian with a De Aquino dynamical operator, and no experimental consequence have been established. The shared integer 5 is exact; its physical interpretation is not.

## 5. Relation to the previous RGS artificial A6/Dirac audit

The repository file experiments/misc/rgs_supergolden_a6.py is a separate experiment on an artificial six-dimensional periodic grid. It uses a candidate symmetrized Dirac operator in Weitzenbock gauge and a scalar shift by sqrt(5)/2 in its De Aquino factor check. Its header explicitly says that equivalence with a covariant teleparallel operator derived from an action has not been demonstrated.

That experiment does not prove that the A5 torus Laplacian is its Dirac operator, that the twelve shortest-dual modes are protected De Aquino zero modes, or that sqrt(5)/2 is forced by A5 geometry. Likewise, the earlier RGS/PSU audit found a nearest Wilson angle about 2.641711, not sqrt(5) about 2.236068 (distance about 0.405643); this document does not identify the Wilson angle with sqrt(5).

## 6. Evidence ledger

| Statement | Status |
|---|---|
| Flat torus, dual metric and Laplacian spectrum | STANDARD / PROVED |
| Minimum k^TQk=5 | PROVED by exact finite certificate |
| First eigenvalue 10*pi^2/3, multiplicity 12 | PROVED |
| Finite graph limit m^2 gamma_m -> 2*pi^2/3 | PROVED in general-modulus analysis |
| De Aquino zero xi^2=5/4 | Algebraic consequence of the stated formula |
| Identification xi_k^2=(k^TQk)/4 | CANDIDATE NORMALIZATION |
| Physical momentum map, operator intertwiner, spectral protection | OPEN |
| Physical validation of the De Aquino expression | NOT ESTABLISHED HERE |

## 7. Reproduction

~~~
PYTHONPATH=. python -m pytest tests/test_gec6_a5_continuum_torus.py -q
PYTHONPATH=. python -m src.gec6_a5_continuum_torus
~~~

Tests cover the exact quadratic form, finite shortest-vector enumeration and bookkeeping of the proposed threshold. They do not establish the physical bridge.


## Revision note: physical normalization developed separately

The earlier formula xi_k^2=(k^T Q k)/4 in this document is only a formal arithmetic match of the shortest-lattice value 5 to De Aquino's threshold 5/4. It is not a physical momentum map and must not be used as one. The dimensionally explicit pentagonal-frame analysis, including the prior balanced variational ansatz, the momentum scale p0, the conditional A5 Compton-length calibration, and the S6 symmetry obstruction, is now documented in [De Aquino–pentagon–A5 momentum normalization](DEAQUINO_PENTAGON_A5_MOMENTUM_NORMALIZATION.md). That newer dossier supersedes the naive matching as the project's current bridge proposal; the physical identification remains open.
