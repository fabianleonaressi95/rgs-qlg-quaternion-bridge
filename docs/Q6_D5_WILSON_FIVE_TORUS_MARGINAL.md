# Q6-D5 Wilson fiber marginal on the annihilator torus

Let C be the integral cycle matrix (129 x 192), and let B be a saturated
integer basis for the rank-124 kernel in cycle coordinates. The edge-flux
fiber lattice is L = B C in Z^192. For beta > 0, set

    f(x) = product_e exp(beta cos(x_e)) / I_0(beta),  x in T^192.

Its Fourier coefficients are the Wilson weights
w(n) = product_e I_|n_e|(beta) / I_0(beta). The annihilator

    H_L = {x in T^192 : ell dot x = 0 mod 2 pi for every ell in L}

therefore gives the exact identity

    W_eff([c]) = integral_H_L f(x) exp(-i (c C) dot x) d Haar(x).

This follows from character orthogonality; absolute convergence of the
Fourier series justifies termwise integration.

For a Smith decomposition U L V = S, with nonzero diagonal entries d_i,
the coordinates x = V z satisfy z_i = 2 pi j_i/d_i for the constrained
coordinates and are free on the remaining 68 coordinates. The finite
component count is product_i d_i. On the Q6 lattice all 124 factors are
one, so H_L is connected and isomorphic to T^68.

The implementation samples only the genuinely free torus coordinates and
any nontrivial finite Smith coordinates. The small scrambled-Sobol run is
a reproducible diagnostic, not a rigorous error interval. A validated
numerical enclosure, a derived pushforward to the D5 Hua measure, and a
physical 3+1D or fine-structure-constant prediction remain open.
