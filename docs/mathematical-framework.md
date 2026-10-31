# Mathematical Framework

## 1. Scope

This document defines the mathematical foundations of the RGS–QLG Quaternion Bridge.

The initial framework concerns:

1. The real quaternion algebra.
2. Unit quaternions and three-dimensional rotations.
3. The relationship between unit quaternions and SU(2).
4. The double covering of SO(3) by SU(2).
5. Composition, inversion and orientation conventions.
6. The extension to discrete lattice holonomies.

The purpose is to establish a consistent mathematical language before introducing broader geometric or physical interpretations.

---

## 2. Quaternion algebra

The quaternion algebra is

\[
\mathbb H =
\{w+x\mathbf i+y\mathbf j+z\mathbf k:
w,x,y,z\in\mathbb R\}.
\]

Its fundamental multiplication rules are

\[
\mathbf i^2=\mathbf j^2=\mathbf k^2=-1,
\]

\[
\mathbf i\mathbf j=\mathbf k,\qquad
\mathbf j\mathbf k=\mathbf i,\qquad
\mathbf k\mathbf i=\mathbf j,
\]

with reversed products changing sign:

\[
\mathbf j\mathbf i=-\mathbf k,\qquad
\mathbf k\mathbf j=-\mathbf i,\qquad
\mathbf i\mathbf k=-\mathbf j.
\]

A quaternion is represented in scalar-first form:

\[
q=(w,x,y,z)=(w,\mathbf v),
\]

where

\[
\mathbf v=(x,y,z).
\]

This component ordering is the canonical convention throughout the repository.

## 3. Quaternion multiplication

Let

\[
q_1=(a,\mathbf u),\qquad
q_2=(b,\mathbf v).
\]

Hamilton multiplication is defined by

\[
q_1q_2=
\left(
ab-\mathbf u\cdot\mathbf v,\,
a\mathbf v+b\mathbf u+\mathbf u\times\mathbf v
\right).
\]

Quaternion multiplication is associative but generally not commutative:

\[
(q_1q_2)q_3=q_1(q_2q_3),
\]

\[
q_1q_2\ne q_2q_1.
\]

All implementations must preserve this multiplication convention.

## 4. Conjugation, norm and inverse

Quaternionic conjugation is

\[
q^*=(w,-\mathbf v).
\]

The squared norm is

\[
|q|^2=qq^*=q^*q
=w^2+x^2+y^2+z^2.
\]

For every nonzero quaternion,

\[
q^{-1}=\frac{q^*}{|q|^2}.
\]

For a unit quaternion,

\[
|q|=1,
\]

and therefore

\[
\boxed{q^{-1}=q^*.}
\]

The implementation must distinguish conjugation from negation:

\[
-q^*=(-w,\mathbf v).
\]

The latter is not the quaternionic inverse.

## 5. Unit quaternions and rotations

The unit quaternions form the group

\[
S^3=
\{q\in\mathbb H:|q|=1\}.
\]

A rotation through angle \(\theta\) around a unit axis

\[
\mathbf n=(n_x,n_y,n_z)
\]

is represented by

\[
q=
\left(
\cos\frac{\theta}{2},
\mathbf n\sin\frac{\theta}{2}
\right).
\]

For a vector \(\mathbf r\in\mathbb R^3\), identify it with the pure quaternion

\[
p=(0,\mathbf r).
\]

The active rotation is

\[
p'=qpq^{-1}.
\]

Its scalar component is zero, and its vector component is the rotated vector.

The factor of one half in the quaternion angle is essential.

## 6. SU(2) representation

Let \(\sigma_x,\sigma_y,\sigma_z\) denote the Pauli matrices.

Associate the unit quaternion

\[
q=(w,x,y,z)
\]

with the complex matrix

\[
U(q)=wI-i(x\sigma_x+y\sigma_y+z\sigma_z).
\]

Explicitly,

\[
U(q)=
\begin{pmatrix}
w-iz & -y-ix\\
y-ix & w+iz
\end{pmatrix}.
\]

For unit quaternions,

\[
U(q)^\dagger U(q)=I,
\]

and

\[
\det U(q)=1.
\]

Thus,

\[
U(q)\in SU(2).
\]

With this convention, quaternion multiplication is represented by matrix multiplication:

\[
U(q_1q_2)=U(q_1)U(q_2).
\]

## 7. The SU(2) to SO(3) double covering

The rotation representation defines a group homomorphism

\[
\rho:SU(2)\longrightarrow SO(3).
\]

Its kernel is

\[
\ker\rho=\{I,-I\}.
\]

Consequently, every three-dimensional rotation has exactly two SU(2) representatives:

\[
\rho(q)=\rho(-q).
\]

The groups are related by

\[
SO(3)\cong SU(2)/\{\pm I\}.
\]

This is a mathematical identity, not a numerical approximation.

## 8. Inverse rotation and antipodal lift

For a unit quaternion,

\[
q^{-1}=q^*.
\]

The two representatives of the inverse rotation are

\[
q^*
\quad\text{and}\quad
-q^*.
\]

They have the same image under the rotation map:

\[
\rho(q^*)=\rho(-q^*).
\]

However, they are distinct elements of SU(2):

\[
qq^*=+1,
\]

\[
q(-q^*)=-1.
\]

This distinction must be preserved in all calculations involving SU(2) composition and lattice transport.

## 9. Composition convention

The active rotation associated with

\[
q=q_1q_2
\]

acts on a pure quaternion as

\[
p'=q_1(q_2pq_2^{-1})q_1^{-1}.
\]

Therefore, the rightmost quaternion acts first.

The corresponding rotation matrices satisfy

\[
R(q_1q_2)=R(q_1)R(q_2).
\]

All source code and tests must use this convention consistently.

## 10. Discrete lattice connections

Consider an oriented graph

\[
G=(V,E).
\]

Assign to each oriented edge \(g\rightarrow h\) a group element

\[
U_{gh}\in SU(2).
\]

The reverse edge is assigned

\[
U_{hg}=U_{gh}^{-1}=U_{gh}^{\dagger}.
\]

For an oriented closed path

\[
\gamma=(v_0,v_1,\ldots,v_n=v_0),
\]

define the holonomy

\[
H_\gamma=
U_{v_0v_1}
U_{v_1v_2}
\cdots
U_{v_{n-1}v_n}.
\]

The multiplication order follows the path convention defined above.

Under a local gauge transformation

\[
g_v\in SU(2),
\]

the link transforms as

\[
U_{vh}\longmapsto
g_v U_{vh}g_h^{-1}.
\]

The holonomy based at \(v_0\) transforms by conjugation:

\[
H_\gamma\longmapsto
g_{v_0}H_\gamma g_{v_0}^{-1}.
\]

Its conjugacy-invariant quantities are therefore gauge invariant.

## 11. Scientific status

The following statements are established mathematical identities:

- Quaternion multiplication rules.
- Quaternionic conjugation and inversion.
- The unit-quaternion representation of rotations.
- The SU(2) matrix representation.
- The SU(2)/SO(3) double covering.
- The composition and inverse relations.
- The transformation law for lattice holonomies.

Their implementation in this repository remains subject to independent numerical verification.

Applications to RGS–PSU, QLG, biverse structures or physical gauge fields require additional mathematical and physical arguments. They are not implied by the identities above.

## 12. Next steps

The implementation will proceed through four modules:

- `quaternions.py`: algebra and unit-quaternion operations.
- `rotations.py`: conversion between quaternions and SO(3) matrices.
- `su2_lifts.py`: SU(2) representation and lift handling.
- `holonomies.py`: discrete parallel transport and closed-loop holonomies.

Each module will be accompanied by independent tests.
