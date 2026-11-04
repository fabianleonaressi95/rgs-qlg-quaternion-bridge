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
