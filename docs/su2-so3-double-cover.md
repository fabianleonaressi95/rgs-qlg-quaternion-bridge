# The SU(2)–SO(3) Double Cover

## 1. Introduction

The groups \(SU(2)\) and \(SO(3)\) describe closely related aspects of three-dimensional rotations.

The group \(SO(3)\) represents proper rotations of Euclidean three-dimensional space.

The group \(SU(2)\) provides a double-covering representation of these rotations.

This relationship is fundamental to quaternionic rotation theory, spinorial representations and the construction of discrete lattice holonomies.

The purpose of this document is to establish the precise mathematical correspondence.

## 2. Definition of SO(3)

The special orthogonal group in three dimensions is

\[
SO(3)=
\{R\in\mathbb R^{3\times3}:
R^TR=I,\ \det R=1\}.
\]

Its elements are real orthogonal matrices representing orientation-preserving rotations.

The group operation is ordinary matrix multiplication.

## 3. Definition of SU(2)

The special unitary group is

\[
SU(2)=
\{U\in\mathbb C^{2\times2}:
U^\dagger U=I,\ \det U=1\}.
\]

Every element can be expressed as

\[
U=
\begin{pmatrix}
a&b\\
-b^*&a^*
\end{pmatrix},
\qquad
|a|^2+|b|^2=1.
\]

Thus, \(SU(2)\) is topologically equivalent to the unit three-sphere \(S^3\).

## 4. Quaternionic representation of SU(2)

Using the canonical scalar-first quaternion

\[
q=(w,x,y,z),
\]

define

\[
U(q)=wI-i(x\sigma_x+y\sigma_y+z\sigma_z).
\]

In matrix form,

\[
U(q)=
\begin{pmatrix}
w-iz&-y-ix\\
y-ix&w+iz
\end{pmatrix}.
\]

For

\[
w^2+x^2+y^2+z^2=1,
\]

we have

\[
U(q)^\dagger U(q)=I
\]

and

\[
\det U(q)=1.
\]

The correspondence preserves multiplication:

\[
U(q_1q_2)=U(q_1)U(q_2).
\]

## 5. Construction of the rotation map

Let \(\mathbf r\in\mathbb R^3\).

Associate to it the traceless Hermitian matrix

\[
X(\mathbf r)=
r_x\sigma_x+r_y\sigma_y+r_z\sigma_z.
\]

For \(U\in SU(2)\), define \(R(U)\) through

\[
X(R(U)\mathbf r)
=
U X(\mathbf r) U^\dagger.
\]

Because unitary conjugation preserves Hermiticity and trace, this transformation defines a real linear map on \(\mathbb R^3\).

It preserves the Euclidean norm, and its determinant is \(+1\).

Consequently,

\[
R(U)\in SO(3).
\]

The map

\[
\rho:SU(2)\rightarrow SO(3)
\]

is a group homomorphism:

\[
\rho(U_1U_2)=\rho(U_1)\rho(U_2).
\]

For the quaternion convention adopted in this repository, this map agrees with the active rotation

\[
p'=qpq^{-1}.
\]

## 6. Axis-angle correspondence

For a unit vector \(\mathbf n\) and angle \(\theta\), define

\[
q=
\left(
\cos\frac{\theta}{2},
\mathbf n\sin\frac{\theta}{2}
\right).
\]

The associated SU(2) element is

\[
U=
\exp\left(
-\frac{i\theta}{2}
\mathbf n\cdot\boldsymbol\sigma
\right).
\]

Its SO(3) image is a rotation through angle \(\theta\) around \(\mathbf n\), using the right-handed active convention.

The half-angle is the essential feature of the double-covering relationship.

## 7. The kernel of the covering map

The kernel is the set of elements mapped to the identity rotation:

\[
\ker\rho=\{I,-I\}.
\]

In quaternionic notation, these elements are

\[
q_+=(1,0,0,0),
\]

\[
q_-=(-1,0,0,0).
\]

Both induce the identity rotation in \(SO(3)\).

They remain distinct elements of \(SU(2)\).

## 8. The two lifts

For every \(R\in SO(3)\), there are exactly two elements \(U,-U\in SU(2)\) such that

\[
\rho(U)=\rho(-U)=R.
\]

In quaternionic notation,

\[
R(q)=R(-q).
\]

The two quaternionic representatives are antipodal points of the unit three-sphere.

The group relation is

\[
SO(3)\cong SU(2)/\{\pm I\}.
\]

This is an exact mathematical identity.

## 9. Quaternionic inversion and the second lift

For a unit quaternion,

\[
q^{-1}=q^*.
\]

The inverse rotation has two SU(2) representatives:

\[
q^*
\qquad\text{and}\qquad
-q^*.
\]

Their images in SO(3) coincide:

\[
R(q^*)=R(-q^*).
\]

However,

\[
qq^*=+1,
\]

whereas

\[
q(-q^*)=-1.
\]

Therefore, the second representative is not the group inverse.

This distinction matters whenever calculations retain information about SU(2), rather than only the corresponding SO(3) rotation.

## 10. Sign consistency in products

Let

\[
q_1,\ldots,q_n\in SU(2).
\]

Replacing individual representatives by their antipodes,

\[
q_j\longmapsto s_jq_j,
\qquad s_j\in\{+1,-1\},
\]

changes the product to

\[
q'_1q'_2\cdots q'_n
=
\left(\prod_{j=1}^{n}s_j\right)
q_1q_2\cdots q_n.
\]

The resulting product may differ by the central element \(-I\).

Nevertheless, both products have the same image under \(\rho\).

Thus, SO(3) rotation data alone cannot determine the global SU(2) sign of a product.

A computational implementation must explicitly record its lift convention.

## 11. Spinorial periodicity

For a rotation about a fixed axis,

\[
U(\theta)=
\exp\left(
-\frac{i\theta}{2}
\mathbf n\cdot\boldsymbol\sigma
\right).
\]

Then

\[
U(2\pi)=-I,
\]

and

\[
U(4\pi)=I.
\]

By contrast, the corresponding SO(3) rotations satisfy

\[
R(2\pi)=I.
\]

This is the standard spinorial \(4\pi\) periodicity of the fundamental representation of \(SU(2)\).

It is a representation-theoretic property, not a claim that ordinary spatial vectors require a \(4\pi\) rotation to return to their original orientation.

## 12. Numerical validation requirements

The implementation must independently verify:

1. \(U^\dagger U=I\).
2. \(\det U=1\).
3. \(\rho(U)\in SO(3)\).
4. \(\rho(U_1U_2)=\rho(U_1)\rho(U_2)\).
5. \(\rho(U)=\rho(-U)\).
6. The kernel consists of \(I\) and \(-I\).
7. \(U(q^*)=U(q)^{-1}\).
8. \(U(-q^*)=-U(q)^{-1}\).
9. The \(2\pi\) and \(4\pi\) periodicity relations.

All tests must report numerical residuals.

## 13. Scope and limitations

The SU(2)–SO(3) double covering is an established mathematical result.

Its use in a discrete lattice model does not automatically establish a physical gauge theory, a particular particle representation or a relationship with the Standard Model.

Those applications require additional assumptions and independent derivations.

The purpose of this repository is first to establish the algebraic and computational foundations with transparent, reproducible tests.

## 14. References

- J. J. Sakurai and J. Napolitano, *Modern Quantum Mechanics*.
- M. Nakahara, *Geometry, Topology and Physics*.
- B. C. Hall, *Lie Groups, Lie Algebras, and Representations: An Elementary Introduction*.
- J. Stillwell, *Naive Lie Theory*.
