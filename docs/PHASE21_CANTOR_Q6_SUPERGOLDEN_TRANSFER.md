# Phase 21 — Cantor–Q6–supergolden coupled transfer operator

Continuation of the canonical 27-chapter RGS–PSU–PRISM–QLG monograph. This is not a new chapter. It formalizes the Cantor/hexagram/supergolden generative core as a finite-state symbolic transfer model.

## Cantor coding and Q6

An infinite binary sequence omega=(b_k) is mapped to the middle-third Cantor set by Pi(omega)=sum_k 2 b_k/3^k. The 64 length-six prefixes are the 64 Q6 vertices. Hamming-distance-one changes form the 6-cube with 64 vertices, 192 edges and degree 6. This identifies finite symbolic prefixes with hexagrams; it does not by itself prove a physical fractal model.

## Supergolden Markov presentation

A_psi = [[1,0,1],[1,0,0],[0,1,0]], with psi^3-psi^2-1=0. The right Perron vector is r=(1,psi^-1,psi^-2). The Parry transition matrix is P_ij=A_ij r_j/(psi r_i), and is row stochastic.

Set q_phi=phi^-1 exp(-phi^-2), g_psi=(r_1 r_2 r_3)^(1/3), q_i=q_phi r_i/g_psi, p_0(i)=1/(1+q_i), p_1(i)=q_i/(1+q_i). This is an explicit model choice coupling Cantor branch weights to the supergolden state.

## Coupled transfer and pressure

(L f)_i(x)=sum_j P_ij sum_{b=0,1} p_b(i) f_j(bx). It preserves constants. A finite state path and emitted bit word have probability pi_i0 product_k P_ik,ik+1 p_b(k+1)(i_k).

The pressure matrix is M_s(i,j)=P_ij^s [p_0(i)^s+p_1(i)^s] on allowed transitions and zero otherwise. At s=0, M_0=2 A_psi. At s=1, M_1=P and its spectral radius is 1.

For the declared synchronous symbolic metric ell_n=(3 psi)^(-n), define D_s=log rho(M_s)/((1-s)log(3 psi)) for s != 1, and D_1=-P'(1)/log(3 psi), where P(s)=log rho(M_s). This is a generalized dimension for the symbolic model and chosen metric. A projected geometric dimension needs a separate proof of separation and overlap control.

## Q6 spectral sector

The module constructs the 64x64 Q6 combinatorial Laplacian and checks its spectrum 2j with multiplicities binomial(6,j). The coincidence between the dimension 3*64=192 and the 192 Q6 edges does not identify the two objects.

No intertwiner with the Q6 Laplacian is assumed. The next task is to construct T and audit epsilon=||T L6-L_acc T||_F/max(1,||T L6||_F).

## Evidence status

Exact identities: supergolden polynomial and matrix identity; Perron/Parry construction; row stochasticity; Q6 graph spectrum; pressure definitions.

Numerically evaluated: stationary distribution, pressure spectral radii, generalized dimensions and residuals in the JSON audit.

Model choices: golden state-dependent branch weights and synchronous metric (3 psi)^(-n).

Open: geometric embedding and overlap control; intertwiner to Q6/Krist-All/PRISM; derivation of psi^(-n) from PSU action; operator convergence; physical observables and experimental discriminants.

## Reproduction

~~~bash
PYTHONPATH=. python -m pytest tests/test_cantor_supergolden_transfer.py -q
PYTHONPATH=. python examples/cantor_supergolden_transfer_audit.py
~~~
