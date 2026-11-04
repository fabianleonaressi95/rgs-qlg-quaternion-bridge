# Q6 coupled action: Euler–Lagrange variations

This note derives the first-variation equations for the finite-dimensional
Euclidean action implemented in \`src/unified_discrete_action.py\`. It does not
identify the model with TEGR, a Dirac discretisation, or a continuum theory.

## 1. Definitions and conventions

Let \(G=(V,E)\) be the oriented Q6 edge set and \(F\) its chosen oriented
plaquette list. At each vertex \(v\), let \(e_v\in GL(4,\mathbb R)\),
\(\psi_{\sigma,v}\in\mathbb C^2\), \(\sigma\in\{+,-\}\), and let
\(U_{uv}\in SU(2)\) be the link on the oriented edge \((u,v)\).

Set
\[
\mathcal V_v=|\det e_v|,\qquad
\rho_v=\|\psi_{+,v}\|^2+\|\psi_{-,v}\|^2,\qquad
\Delta_v=\psi_{+,v}-\psi_{-,v},
\]
\[
d_{\sigma,uv}=\psi_{\sigma,v}-U_{uv}\psi_{\sigma,u},\qquad
D_{uv}=\sum_{\sigma=\pm}\|d_{\sigma,uv}\|^2,\qquad
w_{uv}=(\mathcal V_u\mathcal V_v)^{1/2}.
\]

The action is
\[
S=S_e+S_W+S_{\rm hop}+S_{\rm twin}+S_4,
\]
with the five terms and coupling normalisation exactly as in the implementation:
\[
S_e=\frac{\kappa_e}{2}\sum_{(u,v)\in E}\|e_u-e_v\|_F^2,
\quad
S_W=\beta\sum_{p\in F}\left(1-\frac12\operatorname{ReTr}H_p\right),
\]
\[
S_{\rm hop}=\frac{\kappa_\psi}{2}\sum_{(u,v)\in E}w_{uv}D_{uv},
\quad
S_{\rm twin}=\frac{\mu}{2}\sum_v\mathcal V_v\|\Delta_v\|^2,
\quad
S_4=\frac{\lambda}{4}\sum_v\mathcal V_v\rho_v^2.
\]

Here \(H_p\) is the ordered product of links (or inverse links) around the
oriented plaquette. The determinant sign of each \(e_v\) is held fixed during
a smooth variation; the formula for \(|\det e_v|\) is differentiable on each
non-degenerate sign sector.

## 2. Coframe equation

For an arbitrary real matrix variation \(\delta e_v\),
\[
\delta\mathcal V_v
=\mathcal V_v\,\operatorname{Tr}(e_v^{-1}\delta e_v)
=\mathcal V_v\,e_v^{-T}:\delta e_v.
\]

Define the local volume pressure
\[
\alpha_v=
\frac{\kappa_\psi}{4}
\sum_{w:\{v,w\}\in E}
\frac{(\mathcal V_v\mathcal V_w)^{1/2}}{\mathcal V_v}D_{vw}
+\frac{\mu}{2}\|\Delta_v\|^2
+\frac{\lambda}{4}\rho_v^2.
\]
In this expression \(D_{vw}\) is evaluated using the stored link in whichever
orientation represents the undirected edge; its value is unchanged by reversing
the edge and replacing \(U\) by \(U^\dagger\).

The coframe first variation is
\[
\delta_{e_v}S=
\left[
\kappa_e\sum_{w:\{v,w\}\in E}(e_v-e_w)
+\alpha_v\mathcal V_v e_v^{-T}
\right]:\delta e_v.
\]
Thus the unconstrained stationary equation at every vertex is
\[
\boxed{\;
\kappa_e\sum_{w:\{v,w\}\in E}(e_v-e_w)
+\alpha_v\mathcal V_v e_v^{-T}=0.
\;}
\]

This equation is the stationarity condition of the *implemented fixed-frame
regulariser plus determinant-weighted matter terms*. It is not an Einstein
equation, a tetrad equation of TEGR, or a local-Lorentz-covariant gravitational
field equation. In particular, \(S_e\) compares coframes in a globally fixed
frame.

## 3. Matter-sheet equations

Use the real variation convention
\[
\delta S=2\operatorname{Re}\sum_v
\left(g_{+,v}^\dagger\delta\psi_{+,v}
+g_{-,v}^\dagger\delta\psi_{-,v}\right)+\cdots.
\]
For each sheet \(\sigma\), the hopping contribution at vertex \(v\) is
\[
g^{\rm hop}_{\sigma,v}
=\frac{\kappa_\psi}{2}
\left[
\sum_{u:(u,v)\in E}w_{uv}d_{\sigma,uv}
-\sum_{w:(v,w)\in E}w_{vw}U_{vw}^\dagger d_{\sigma,vw}
\right].
\]
The sums are over incoming and outgoing oriented edges in the stored
orientation. The on-site terms give
\[
g^{\rm twin}_{+,v}=\frac{\mu}{2}\mathcal V_v\Delta_v,\qquad
g^{\rm twin}_{-,v}=-\frac{\mu}{2}\mathcal V_v\Delta_v,
\]
\[
g^{(4)}_{\sigma,v}=\frac{\lambda}{2}\mathcal V_v\rho_v\psi_{\sigma,v}.
\]
Consequently
\[
\boxed{\quad
g_{\sigma,v}:=g^{\rm hop}_{\sigma,v}
+g^{\rm twin}_{\sigma,v}
+g^{(4)}_{\sigma,v}=0,\qquad \sigma=\pm.
\quad}
\]
These are complex stationary equations under independent variations of the
two doublets. They are graph-covariant bosonic equations for the chosen
quadratic hopping energy, not fermionic Dirac equations.

## 4. Link equation on \(SU(2)\)

A tangent variation of a stored link can be written
\[
U_{uv}(t)=e^{tX_{uv}}U_{uv},\qquad
X_{uv}^\dagger=-X_{uv},\quad\operatorname{Tr}X_{uv}=0.
\]
Writing \(q_{\sigma,uv}=U_{uv}\psi_{\sigma,u}\), the hopping variation is
\[
\delta_U S_{\rm hop}
=-\kappa_\psi w_{uv}
\operatorname{Re}\sum_{\sigma=\pm}
d_{\sigma,uv}^\dagger X_{uv}q_{\sigma,uv}.
\]

For every plaquette containing this link, differentiate its ordered
holonomy by the product rule. If its ordered factors are \(L_1,\ldots,L_4\),
then
\[
\delta H_p=\sum_{j=1}^4
L_1\cdots L_{j-1}(\delta L_j)L_{j+1}\cdots L_4,
\]
where \(\delta L_j=X_{uv}U_{uv}\) for a forward occurrence and
\(\delta L_j=-U_{uv}^\dagger X_{uv}\) for an inverse occurrence.
The Wilson variation is
\[
\delta_U S_W=-\frac{\beta}{2}
\sum_{p\ni (u,v)}\operatorname{ReTr}(\delta H_p).
\]
Therefore the exact link stationarity condition is the weak equation
\[
\boxed{\quad
-\frac{\beta}{2}\sum_{p\ni (u,v)}
\operatorname{ReTr}(\delta_X H_p)
-\kappa_\psi w_{uv}\operatorname{Re}
\sum_{\sigma=\pm}d_{\sigma,uv}^\dagger X_{uv}q_{\sigma,uv}=0
\quad\text{for every }X_{uv}\in\mathfrak{su}(2).
\quad}
\]
This is the vanishing of the left-trivialised \(SU(2)\) gradient. It avoids
assuming a coordinate chart for the group and includes both forward and
backward occurrences in plaquette paths.

## 5. Required numerical verification before adding a solver

An implementation of these gradients must pass independent directional
finite-difference tests, not merely compare two code paths sharing the same
formula:

1. Coframe: compare the analytic contraction with
   \([S(e+\epsilon h)-S(e-\epsilon h)]/(2\epsilon)\), keeping all determinants
   nonzero and their signs fixed.
2. Matter: compare \(2\operatorname{Re}(g^\dagger h)\) with the same centred
   difference for independent complex perturbations of each sheet.
3. Link: compare the weak \(\mathfrak{su}(2)\) derivative with a centred
   difference along \(e^{\pm\epsilon X}U\), including links that occur in a
   plaquette in both orientations.
4. Sweep \(\epsilon\) over several decades and report the truncation/roundoff
   curve, rather than one tolerance at one step size.
5. Recheck local \(SU(2)\) gauge invariance after any gradient or solver
   implementation.

No Euler–Lagrange solver, gradient code, or finite-difference validation is
claimed by this derivation note. The next implementation gate is to encode
these derivatives and test them independently.

## 6. Boundary of the present model

The present link \(U_{uv}\) is an internal \(SU(2)\) gauge transporter. It is
not derived from the coframe. No spin connection, Clifford module, Lorentzian
signature, torsion scalar, or TEGR boundary term appears in this action.
A future gravity/fermion extension must introduce those structures explicitly
and prove their transformation laws; they cannot be inferred from the
stationarity equations above.
