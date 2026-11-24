# De Aquino momentum normalization: pentagonal frame and the A5 scale

**Purpose.** Recover the De Aquino momentum variable from the earlier source audit, put the regular-pentagon frame on an explicit variational footing, and identify exactly what is and is not derived when comparing it with the A5 continuum torus.

**Evidence labels:** algebraic consequence / geometric theorem / variational ansatz / physical identification open.

## 1. De Aquino variable: keep its physical dimensions

The formula transcribed from De Aquino (2002) in the earlier project audit is
\[
\frac{m_g}{m_{i0}}=3-2\sqrt{1+\xi^2},\qquad
\xi=\frac{\Delta p}{m_{i0}c}
=\frac{n_r U}{m_{i0}c^2},
\]
where the same source convention writes \(\Delta p=n_rU/c\). Here \(n_rU\) has dimensions of energy, \(\Delta p\) of momentum, and \(\xi\) is dimensionless. The zero of the displayed mass-ratio expression is
\[
3=2\sqrt{1+\xi^2}
\quad\Longrightarrow\quad
\boxed{\xi_*^2=\frac54,\quad |\xi_*|=\frac{\sqrt5}{2}.}
\]
This is an exact algebraic threshold of the stated formula, not independent experimental validation of that formula.

## 2. The regular-pentagon frame

Let
\[
e_j=(\cos(2\pi j/5),\sin(2\pi j/5)),\quad
E=(e_0\ e_1\ e_2\ e_3\ e_4)\in\mathbb R^{2\times5}.
\]
The tight-frame identity is
\[
EE^T=\frac52 I_2.
\]
Consequently
\[
P:=\frac25E^TE
\]
is an orthogonal projector of rank two on \(\mathbb R^5\), and
\[
G_5:=E^TE=\frac52P,\qquad
\operatorname{spec}(G_5)=\left(\frac52,\frac52,0,0,0\right).
\]
For a unit five-channel amplitude \(q\in S^4\), write
\[
q=q_{\rm att}+q_{\rm null},\quad
q_{\rm att}=Pq,\quad q_{\rm null}=(I-P)q.
\]
Then
\[
q^TG_5q=\|Eq\|^2=\frac52\|q_{\rm att}\|^2.
\]
The De Aquino threshold is equivalent to
\[
q^TG_5q=\frac54
\quad\Longleftrightarrow\quad
\|q_{\rm att}\|^2=\|q_{\rm null}\|^2=\frac12.
\]
The threshold set is \(S^1(1/\sqrt2)\times S^2(1/\sqrt2)\subset S^4\). The frame identity alone does not select this set or an individual \(q\).

## 3. Balanced split from the previously proposed variational ansatz

The prior discussion proposed, rather than derived from De Aquino dynamics, the constrained functional
\[
\mathcal E(q)=-\lambda\|Pq\|^2\|(I-P)q\|^2,
\qquad \lambda>0,\quad \|q\|=1.
\]
Put \(x=\|Pq\|^2\). Orthogonality gives \(\|(I-P)q\|^2=1-x\), so
\[
\mathcal E(q)=-\lambda x(1-x)
=\lambda\left(x-\frac12\right)^2-\frac\lambda4.
\]
Thus the global minimum of this chosen functional occurs at \(x=1/2\). It follows algebraically that
\[
q^TG_5q=\frac52x=\frac54.
\]
An equivalent auxiliary-field ansatz is
\[
F(a,b,Y)=F_0+\frac\kappa2\|Y\|^2-g\langle Y,a\otimes b\rangle,
\quad a=Pq,\quad b=(I-P)q,\quad \kappa>0.
\]
Eliminating the stationary auxiliary field \(Y_*=g(a\otimes b)/\kappa\) gives
\[
F_{\rm eff}=F_0-\frac{g^2}{2\kappa}\|a\|^2\|b\|^2,
\]
the same form with \(\lambda=g^2/(2\kappa)>0\).

**Interpretation limit:** this proves the balanced split is the minimizer of that auxiliary-field model. It does not derive the coupling, auxiliary field, or functional from De Aquino's source equation, the A5 Laplacian, or an independently established action. It is a candidate variational mechanism, not yet a physical theorem.

## 4. Restore the physical momentum scale

A unit vector \(q\) is dimensionless. The frame output \(Eq\) is also dimensionless. Introduce a momentum scale explicitly:
\[
\boxed{\Delta p_{\rm pent}:=p_0\|Eq\|,\qquad
\xi_{\rm pent}:=\frac{\Delta p_{\rm pent}}{m_{i0}c}
=\frac{p_0}{m_{i0}c}\sqrt{q^TG_5q}.}
\]
Here \(p_0\) is a physical momentum scale to be fixed by the model. At the balanced variational minimum,
\[
\xi_{\rm pent}^2
=\left(\frac{p_0}{m_{i0}c}\right)^2\frac54.
\]
Therefore the exact De Aquino threshold is recovered if and only if
\[
\boxed{p_0=m_{i0}c.}
\]
This is the rest-momentum scale (equivalently \(\hbar/\bar\lambda_C\), with reduced Compton wavelength \(\bar\lambda_C=\hbar/(m_{i0}c)\)). It is a dimensionally natural calibration, but the pentagon frame and variational ansatz do not independently force it. It must be derived from a physical coupling or fixed by an independently specified boundary/dispersion condition.

Using the source convention \(\Delta p=n_rU/c\), the same calibration reads \(n_rU=m_{i0}c^2\). Without imposing this rest-momentum calibration, the De Aquino threshold itself requires \(n_rU=m_{i0}c^2\sqrt5/2\). Do not conflate these two equations.

## 5. Independent A5 continuum momentum and its required scale

For the flat torus \(T_{A_5}=V/A_5\), the dual mode labelled by \(k\in\mathbb Z^5\) has dimensionless dual norm
\[
\|\xi_k\|^2=\frac16 k^TQk,\qquad Q=6C_{A_5}^{-1}.
\]
If the torus is assigned a physical length scale \(L\), its Fourier momentum magnitude is
\[
p_k=\frac{2\pi\hbar}{L}\|\xi_k\|
=\frac{2\pi\hbar}{L}\sqrt{\frac{k^TQk}{6}}.
\]
For the shortest nonzero modes, \(k^TQk=5\), hence
\[
\frac{p_{\min}}{m_{i0}c}
=\frac{2\pi\bar\lambda_C}{L}\sqrt{\frac56}.
\]
Demanding that this A5 mode itself equal the De Aquino threshold gives the necessary length calibration
\[
\frac{2\pi\bar\lambda_C}{L}\sqrt{\frac56}
=\frac{\sqrt5}{2}
\quad\Longleftrightarrow\quad
\boxed{L=\frac{4\pi}{\sqrt6}\bar\lambda_C.}
\]
This is a conditional calibration equation, not a prediction of the unscaled mathematical torus: the lattice quotient fixes shape and dimensionless spectrum, not physical length. A mechanism fixing \(L\) to this Compton-scale value remains to be supplied.

## 6. Missing map between A5 modes and pentagon amplitudes

The two constructions now have explicit dimensionally correct formulas, but are not yet the same construction:
- A5 supplies a character label \(k\), a dual norm \(k^TQk/6\), and a Fourier momentum after a physical length \(L\) is specified.
- The pentagon supplies a dimensionless five-channel amplitude \(q\), a rank-two projector \(P\), and the projected magnitude \(\|Eq\|\).
- A bridge requires a specified map \(J:A_5^*\to S^4\) (or a corresponding map on mode spaces) and a physical identification of \(p_k\) with \(p_0\|EJ(k)\|\).

There is a symmetry obstruction to a **canonical fully S6-equivariant** identification of the pentagon plane with an A5 subspace. The A5 standard five-dimensional representation of S6 is irreducible. Any operator commuting with all of S6 on this real irreducible representation is scalar (Schur's lemma); a rank-two projector \(P\) cannot commute with the full S6 action. The regular pentagon frame has dihedral symmetry \(D_5\), not the full S6 symmetry of the A5 root geometry. A chosen pentagon embedding therefore introduces extra structure and breaks full S6 symmetry. A map may still be built after choosing that structure, but it cannot be called canonical and S6-equivariant.

The earlier formal expression \(\xi_k^2=(k^TQk)/4\) merely matches threshold numbers. It does not supply the missing map or physical scale and is superseded here by the dimensionally explicit equations above.

## 7. Evidence ledger and next theorem target

| Statement | Status |
|---|---|
| De Aquino threshold \(\xi_*^2=5/4\) from the quoted formula | ALGEBRAIC CONSEQUENCE |
| Pentagon tight-frame identity, rank-two projector and spectrum | EXACT GEOMETRY |
| Balanced split minimizes the stated \(\mathcal E(q)\) | PROVED FOR THE ANSATZ |
| Auxiliary-field elimination yields the product functional | ALGEBRAIC, GIVEN THE ANSATZ |
| Momentum conversion \(\Delta p=p_0\|Eq\|\) | DEFINITION OF CANDIDATE PHYSICAL MAP |
| Required calibration \(p_0=m_{i0}c\) | CONDITIONAL NECESSARY NORMALIZATION |
| A5 physical length calibration \(L=4\pi\bar\lambda_C/\sqrt6\) | CONDITIONAL CONSEQUENCE |
| Canonical map \(A_5^*\to S^4\) compatible with the pentagon | OPEN; full S6-equivariance obstructed |
| Derivation of variational coupling and physical scale from a common action | OPEN |
| Independent physical validation of De Aquino's mass-ratio law | NOT ESTABLISHED |

**Next theorem target:** construct an explicitly symmetry-broken, physically dimensioned map \(J\) from A5 dual modes to five pentagonal channel amplitudes, derive its scale and coupling from a common action (rather than choosing them to fit \(\sqrt5/2\)), and test whether the resulting operator reproduces both the A5 continuum eigenvalue and De Aquino's momentum dependence. Until then, the bridge remains a structured candidate, not a proved physical identification.
