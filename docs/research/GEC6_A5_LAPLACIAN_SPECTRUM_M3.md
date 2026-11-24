# GEC-6 spectral audit: the finite A5-root Cayley graph at m=3

**Status:** exact Fourier spectrum, independently confirmed by numerical diagonalization of the full 243-by-243 Laplacian. No physical or continuum equivalence is claimed.

## 1. Graph and Laplacian

Let
\[
G=(\mathbb Z/3\mathbb Z)^5,\qquad
R=\{e_i-e_j:i\ne j\}\subset A_5/3A_5.
\]
The graph \(\Gamma_3=\operatorname{Cay}(G,R)\) is connected, simple, undirected and 30-regular, with 243 vertices and 3,645 edges. Its combinatorial Laplacian is
\[
L=30I-A.
\]

## 2. Exact Fourier derivation

The character group of \(G\) is again \((\mathbb Z/3\mathbb Z)^5\). In simple-root coordinates, a character indexed by \(k\) evaluates a root step \(r\) as
\[
\chi_k(r)=\exp\!\left(\frac{2\pi i}{3}\langle k,r\rangle\right).
\]
Pair the 30 roots into the 15 opposite pairs \(\{r,-r\}\). For each pair, the adjacency contribution is
\[
\chi_k(r)+\chi_k(-r)=
\begin{cases}
2,&\langle k,r\rangle=0\pmod3,\\
-1,&\langle k,r\rangle\ne0\pmod3.
\end{cases}
\]
Let \(a(k)\) be the number of the 15 positive-root pairs with zero phase. Then
\[
\lambda_A(k)=2a(k)-(15-a(k))=3a(k)-15,
\]
and therefore
\[
\boxed{\lambda_L(k)=30-\lambda_A(k)=45-3a(k).}
\]
This is an exact integer formula; no floating-point trigonometry is needed for the Fourier spectrum.

Enumerating all 243 characters gives:

| Laplacian eigenvalue | Multiplicity |
|---:|---:|
| 0 | 1 |
| 15 | 12 |
| 24 | 30 |
| 27 | 50 |
| 33 | 120 |
| 36 | 30 |
| **Total** | **243** |

Thus
\[
\boxed{\operatorname{Spec}(L_{\Gamma_3})
=\{0^{(1)},15^{(12)},24^{(30)},27^{(50)},33^{(120)},36^{(30)}\}.}
\]
The graph is connected, so the zero eigenvalue is simple. The (unnormalized) spectral gap is 15. Since the graph is 30-regular, the normalized Laplacian is \(L/30\), with gap
\[
\boxed{\lambda_1(L/30)=\tfrac12.}
\]

## 3. Independent matrix check

The audit runner constructs the complete symmetric integer adjacency matrix \(A\in\{0,1\}^{243\times243}\), checks \(A=A^\top\) and row sums 30, and diagonalizes \(30I-A\) using a symmetric eigensolver. All 243 numerical eigenvalues lie within \(10^{-9}\) of integers, and the rounded multiplicities agree exactly with the 243-character Fourier count.

This numerical diagonalization is an independent implementation check. The exact multiplicities themselves follow from the finite abelian Fourier decomposition and the integer root-pair formula above.

## 4. Comparison with Q6 and PRISM

Compare normalized Laplacians to remove the trivial degree scaling.

| Graph | Degree | Normalized spectral gap |
|---|---:|---:|
| GEC-6 \(\Gamma_3\) | 30 | \(0.5\) |
| Q6 hypercube \(Q_6\) | 6 | \(1/3\) |
| PRISM \(C_{64}\square C_6^3\square C_5\) | 10 | \((2-2\cos(2\pi/64))/10\approx0.000963054666\) |

For Q6, the combinatorial Laplacian spectrum is \(2k\) with multiplicity \(\binom6k\), \(k=0,\ldots,6\). Its normalized gap is \(2/6=1/3\).

For the Cartesian PRISM product, the Laplacian spectrum is the sum of the five cycle Laplacian eigenvalues. The smallest nonzero eigenvalue comes from the \(C_{64}\) factor:
\[
\lambda_1(P)=2-2\cos(2\pi/64)\approx0.009630546656.
\]
The degree is 10, giving the normalized gap shown in the table.

The finite GEC-6 quotient has a normalized gap larger than both Q6 and the PRISM product at this particular modulus. This is a comparison of three finite graphs, not a claim that the family is an expander family.

## 5. Important limit and interpretation

The gap \(1/2\) at \(m=3\) does not mean the family \(\Gamma_m\) has a gap bounded away from zero as \(m\to\infty\). For fixed rank-five lattice generators and larger periodic quotients, low-frequency characters approach the trivial character and the normalized spectral gap scales to zero on the order of \(m^{-2}\). The finite \(m=3\) quotient is a coarse periodic discretization, not evidence of a continuum spectral gap.

The spectrum does not establish:
- graph isomorphism or spectral equivalence with the original PRISM or Q6;
- a canonical map between their eigenvectors or eigenspaces;
- a physical interpretation, field equation, or continuum limit;
- a uniform spectral gap over all moduli \(m\).

## 6. Evidence ledger

| Statement | Status |
|---|---|
| Character formula \(\lambda_L(k)=45-3a(k)\) | PROVED |
| Multiplicity table for \(m=3\) | EXACT FINITE ENUMERATION |
| Agreement with dense 243-by-243 numerical diagonalization | NUMERICAL CROSS-CHECK |
| Normalized gap comparison with Q6 and PRISM | EXACT FORMULAS / NUMERICAL VALUES |
| Uniform expansion as \(m\to\infty\) | FALSE; gap tends to zero |
| Spectral/dynamical/physical equivalence with PRISM | OPEN; not claimed |

## 7. Reproduction

\`\`\`bash
PYTHONPATH=. python -m pytest tests/test_gec6_a5_spectral.py -q
PYTHONPATH=. python examples/gec6_a5_spectral_audit.py
\`\`\`



## 8. General modulus and the asymptotic gap

For every integer m>=3, the normalized Laplacian character eigenvalue is
\[
\mu_m(k)=\frac1{15}\sum_{0\le i<j\le5}
\left(1-\cos\left(\frac{2\pi\langle k,r_{ij}\rangle}{m}\right)\right),
\qquad k\in(\mathbb Z/m\mathbb Z)^5.
\]
The normalized gap is \(\gamma_m=\min_{k\ne0}\mu_m(k)\). The implementation enumerates all \(m^5\) characters, so the finite values below are exhaustive for each listed modulus. Cosines are evaluated numerically; unlike the m=3 table, these decimal values are not claimed to be exact integer spectra.

### Finite-modulus gap audit

| m | Vertices | Normalized gap | m^2 gap |
|---:|---:|---:|---:|
| 3 | 243 | 0.500000000000 | 4.500000000000 |
| 4 | 1,024 | 0.333333333333 | 5.333333333333 |
| 5 | 3,125 | 0.230327668542 | 5.758191713542 |
| 6 | 7,776 | 0.166666666667 | 6.000000000000 |
| 7 | 16,807 | 0.125503399380 | 6.149666569641 |
| 8 | 32,768 | 0.097631072938 | 6.248388668020 |
| 9 | 59,049 | 0.077985185627 | 6.316800035787 |
| 10 | 100,000 | 0.063661001875 | 6.366100187502 |

The sequence \(m^2\gamma_m\) approaches \(2\pi^2/3=6.579736267393\ldots\).

### Proof of the asymptotic constant

In simple-root coordinates, the positive roots \(r_{ij}\) satisfy
\[
Q=\sum_{i<j}r_{ij}r_{ij}^{T}=
\begin{pmatrix}
5&4&3&2&1\\
4&8&6&4&2\\
3&6&9&6&3\\
2&4&6&8&4\\
1&2&3&4&5
\end{pmatrix}.
\]
The minimum of \(k^TQk\) over nonzero integer \(k\) is 5, attained for example at \(k=(1,0,0,0,0)\). Here is a finite exact certificate: the leading principal minors of \(Q-\frac32I\) are \(7/2,27/4,81/8,189/16,243/32\), so \(Q\succ\frac32I\). If an integer vector had \(k^TQk<5\), it would therefore have \(\|k\|^2<10/3\), hence \(\|k\|^2\le3\). Exhaustively checking the finite set of integer vectors with squared norm 1, 2 or 3 gives minimum 5. This establishes the shortest-vector claim without an unbounded search.

For fixed nonzero integer \(k\), Taylor expansion gives
\[
\mu_m(k)=\frac{2\pi^2}{15m^2}k^TQk+O(m^{-4}).
\]
Choosing a shortest vector gives the upper bound. Conversely, any sequence of minimizing characters with \(m^2\mu_m\) bounded must approach the trivial character in each simple-root coordinate (the five simple roots occur among the generators); centered representatives then admit the same Taylor expansion, and the nonzero integer quadratic form is at least 5. Therefore
\[
\boxed{\lim_{m\to\infty}m^2\gamma_m=\frac{2\pi^2}{15}\,5=\frac{2\pi^2}{3}.}
\]

This is a vanishing-gap family, not an expander family. The result concerns only the normalized graph Laplacian of the A5-root Cayley quotients; it adds no spectral or physical equivalence with Q6 or the original PRISM product.
