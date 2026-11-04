# Q6 U(1) character expansion and integer electric flux

## Exact character expansion

For each temporal plaquette angle (F_{0e}), the compact U(1) Boltzmann
factor has the exact Fourier-Bessel expansion

[
e^{\beta_t\cos F_{0e}}=
\sum_{n_e\in\mathbb Z} I_{n_e}(\beta_t)e^{i n_e F_{0e}},
\qquad I_{-n}=I_n.
]

Writing the temporal flux as
(F_{0e}(t)=A_e(t)+A_0(b,t)-A_e(t+1)-A_0(a,t)), the coefficient of
(A_0(v,t)) in the Fourier phase is the oriented incidence divergence
((D_0^T n(t))_v). Normalized Haar integration therefore gives

[
\int_0^{2\pi}\frac{dA_0(v,t)}{2\pi}
 e^{i A_0(v,t)(D_0^Tn(t))_v}
=\delta_{(D_0^Tn(t))_v,0}.
]

Thus the exact electric-flux constraint is

[
n(t)\in\mathbb Z^{192},\qquad D_0^Tn(t)=0.
]

## Q6 cycle lattice

The Q6 graph is connected, with 64 vertices and 192 oriented edges.
The incidence matrix has rank 63, so its integer cycle lattice has rank
(192-63=129). The implementation constructs a fundamental integral
cycle basis from a spanning tree and verifies the divergence identity with
integer arithmetic (not a floating-point tolerance).

The temporal character weight for a link carrying integer flux (n_e) is
(I_{n_e}(\beta_t)). The module computes a controlled finite Fourier
truncation for numerical audits; the exact expansion itself is untruncated.

## Scope boundary

This is an exact finite-lattice U(1) dual electric-flux representation on
(Q6\times C_{N_t}), which is 6+1-dimensional Euclidean. It does not
derive a physical 3+1D reduction, a continuum limit, the anisotropic
scaling of (\beta_s,\beta_t), or their derivation from the RGS-PSU action.
In particular, the Maxwell-normalized five-mode Hamiltonian remains a
conditional normalization until it is derived from this Wilson transfer
operator with explicit temporal-spacing conventions.
