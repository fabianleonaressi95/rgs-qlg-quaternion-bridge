# Exact annihilator-torus formula for the Q6-D5 Wilson fiber marginal

Let C be the integral cycle matrix (129 x 192), K_Z the saturated rank-124
integer kernel in cycle coordinates, and B_K a full integral row basis for
K_Z. Its image in edge-flux coordinates is the rank-124 row lattice
L = B_K C subset Z^192.

For beta>0 define
f(x) = product_e exp(beta cos(x_e))/I_0(beta), x in T^192.
Its Fourier coefficients are
fhat(n) = product_e I_|n_e|(beta)/I_0(beta) = w(n).

Let H_L be the annihilator subgroup
H_L = {x in T^192 : ell dot x = 0 mod 2pi for all ell in L}.
Character orthogonality on this compact subgroup gives the exact identity
W_eff([c]) = integral_{H_L} f(x) exp(-i(cC).x) dmu_Haar(x).
Indeed, the integral of exp(i(n-cC).x) is one exactly when n-cC belongs
to L, and zero otherwise. Since f has an absolutely convergent Fourier
series, termwise integration is justified. This is the complete fiber sum.

The annihilator has dimension 192-124=68. It need not be connected.
If U L V = S is a Smith decomposition with V unimodular and nonzero
diagonal invariants d_1,...,d_124, write x=V z. Then
z_i=2pi*j_i/d_i for i<124, j_i=0,...,d_i-1, while the last 68
coordinates are free torus variables. Haar measure is the uniform average
over all finite components and normalized Haar measure on T^68.

The implementation computes this parameterization exactly using Smith
normal form. Its optional scrambled-Sobol evaluation is explicitly a QMC
estimate, not a certified interval. The exact Fourier identity is proved;
numerical certification and a Wilson-to-Hua measure identification remain
open. The representation is not a five-dimensional integral: that tempting
reduction would incorrectly identify w(cC) with the Fourier coefficients
of a composed product-cosine function.


## Numerical parameterization refinement (PR #160)

The exact Smith audit on the merged Q6 lattice gives 124 invariant factors,
all equal to one. Thus the embedded fiber lattice L is primitive in Z^192:
the annihilator is connected, with no finite component factor, and is
isomorphic to T^68.

The SNF right transform V gives x=V z. Since all Smith diagonal entries
are one, the first 124 coordinates of z vanish and only the final 68 are
free. The QMC implementation now samples these 68 free coordinates directly.
The previous implementation used a 192-dimensional Sobol sequence while
124 coordinates were deterministic zeros; that wasted low-discrepancy
coordinates and produced a highly unstable pilot estimate (relative
between-scramble standard error about 0.73 at 256 samples/scramble).

The revised CI diagnostic uses 2^12 points per scramble and eight independent
scrambles. It reports the intrinsic dimension and every scramble estimate.
This is a convergence diagnostic, not a validated error interval.

The D5 Hua density remains a distinct object. A quantitative comparison
requires a specified map from the five quotient charges to the D5 generic
norm and a derived pushforward measure. No such map is asserted here.
