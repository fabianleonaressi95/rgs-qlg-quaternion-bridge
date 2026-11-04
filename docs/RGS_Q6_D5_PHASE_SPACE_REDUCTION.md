# Q6 spectral phase space and the type-IV_5 tangent representation

Status: exact finite-dimensional spectral construction plus numerical
representation audits. This is the first non-arbitrary 10-real-dimensional
candidate sector. It is **not** a global identification of Q6 with Wyler's
bounded domain and it does not derive alpha.

## 1. Start from the gauge quotient, not from ten selected coordinates

The compact Q6 U(1) link quotient has 129 real configuration coordinates.
Its tangent at the identity is the Coulomb-gauge space
[
\mathcal H^1_{\rm coex}=\ker D_0^T\subset C^1(Q6),\qquad
\dim\mathcal H^1_{\rm coex}=192-63=129.
]
Use the cubical Hodge 1-Laplacian
[
\Delta_1=D_0D_0^T+D_1^TD_1.
]
On the quotient tangent, the exact part is removed and the coexact operator
is (D_1^TD_1). Its spectrum is
[
4^{[15]},\quad6^{[40]},\quad8^{[45]},\quad10^{[24]},\quad12^{[5]}.
]
The top eigenvalue 12 is isolated and has multiplicity exactly five. Its
spectral projector is independent of any eigenvector basis. The implementation
extracts a deterministic orthonormal frame from that projector.

## 2. The selected five-dimensional space is the A5 irrep 5

The 60 proper icosahedral rotations act on Q6 cochains and commute with the
cubical Hodge Laplacian. Restricting the induced action to the top eigenspace
gives the A5 character, ordered by conjugacy classes
[
(1,15,20,12,12):qquad (5,1,-1,0,0).
]
Its inner products with the five irreducible A5 characters give multiplicities
[
m_1=m_3=m_{3'}=m_4=0,qquad m_5=1.
]
Thus this is one irreducible real 5-dimensional representation, not a
10-dimensional subspace chosen by hand.

## 3. The ten-dimensional object is its canonical cotangent phase space

Let (V_5) be the selected mode space. Its cotangent lift is
[
T^*V_5\simeq V_5\oplus V_5^*,\qquad \dim_{\mathbb R}=10.
]
Use the inherited Euclidean mode metric to identify (V_5^*\simeq V_5).
In canonical coordinates ((q,p)\in\mathbb R^5\oplus\mathbb R^5),
[
\Omega=\begin{pmatrix}0&I\\-I&0\end{pmatrix}.
]
The quadratic Wilson potential on this eigenspace is
[
V^{(2)}(q)=\frac{\beta\,12}{2}\,q^Tq.
]
The electric kinetic term is taken to be (T(p)=\frac12p^Tp). This unit
kinetic metric is an explicit Hamiltonian convention; it is not yet derived
from a physical Maxwell normalization. The oscillator frequency is
(omega=\sqrt{12\beta}).

The compatible complex structure and positive Kähler metric are
[
J=\begin{pmatrix}0&-I/\omega\\\omega I&0\end{pmatrix},
\qquad
g=\Omega J=\begin{pmatrix}\omega I&0\\0&I/\omega\end{pmatrix}.
]
They obey (J^2=-I), (J^T\Omega J=\Omega), and (g>0).
The phase rotations (e^{\vartheta J}\) give an (SO(2)) symmetry. The
isotropic five-mode oscillator has (SO(5)) symmetry, with the icosahedral
(A_5\subset SO(5)) acting in its irreducible 5 representation. Thus the
linear isotropy representation is
[
\mathbb R^5\otimes\mathbb R^2\simeq\mathbb C^5,
]
with real dimension ten and complex dimension five.

## 4. Match to the tangent of Cartan's type-IV domain

The standard bounded realization is
[
D^{IV}_5=\{z\in\mathbb C^5:
1-2\,z^\dagger z+|z^Tz|^2>0,\quad |z^Tz|<1\}.
]
Its automorphism group has identity component locally
(SO_0(5,2)), with isotropy (SO(5)\times SO(2)). The tangent at the
origin is (mathbb C^5), carrying exactly the (SO(5)\times SO(2))
isotropy representation above.

For the Bergman-kernel convention
[
K(z,z)=C\,[1-2z^\dagger z+|z^Tz|^2]^{-5},
qquad
g_{i\bar j}=\partial_i\partial_{\bar j}\log K,
]
the Hermitian metric at the origin is (g_{i\bar j}(0)=10\delta_{ij}).
The local coordinate normalization is therefore
[
z_D=\frac{1}{\sqrt{10}}\,
\frac{\sqrt\omega\,q+i p/\sqrt\omega}{\sqrt2}.
]
The implementation checks membership in the type-IV domain for local
coordinates and verifies the symplectic, complex, Kähler and isotropy
identities.

## 5. What has and has not been constructed

**Verified finite structure**
- the 129-dimensional Coulomb-gauge tangent;
- the isolated five-dimensional top coexact Hodge eigenspace;
- its irreducible A5 character ((5,1,-1,0,0));
- the 10-dimensional cotangent phase-space lift;
- a compatible complex structure, symplectic form and positive Kähler metric;
- the (SO(5)\times SO(2)) linear oscillator isotropy, with A5 as a concrete
  finite subgroup;
- equality of the resulting tangent representation with that of (D^{IV}_5).

**Still open**
- deriving the unit electric kinetic metric from the RGS/PSU action;
- extending the local tangent identification to a global nonlinear map into
  (D^{IV}_5);
- proving a push-forward measure or Jacobian from the compact Q6 quotient;
- deriving Wyler's finite Euclidean coordinate volume or its radius from RGS;
- physical charge normalization, continuum QED limit and alpha.

The crucial distinction: the compact 129-dimensional gauge quotient is not
diffeomorphic to the noncompact D5 domain. What is constructed is a canonical
10-dimensional **phase-space sector and tangent-representation match**, not
a global quotient equivalence.

## Reproduction
```bash
PYTHONPATH=. python -m pytest tests/test_rgs_q6_d5_phase_space.py -q
PYTHONPATH=. python -c "from src.rgs_q6_d5_phase_space import build_q6_d5_phase_space_sector; print(build_q6_d5_phase_space_sector()['all_checks_passed'])"
```
