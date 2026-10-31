# Quaternion Conventions

## 1. Purpose

This document establishes the canonical quaternion conventions used throughout the RGS–QLG Quaternion Bridge.

All mathematical derivations, Python modules, numerical experiments and tests must follow these conventions unless an alternative convention is explicitly declared.

The primary objective is to guarantee consistency between:

- Quaternion algebra.
- Three-dimensional active rotations.
- SU(2) matrix representations.
- SO(3) rotation matrices.
- Discrete parallel transport.
- Lattice holonomies.

---

## 2. Component ordering

The canonical quaternion representation is scalar-first:

\[
\boxed{q=(w,x,y,z)}
\]

Equivalently,

\[
q=w+x\mathbf i+y\mathbf j+z\mathbf k.
\]

The scalar part is

\[
\operatorname{Re}(q)=w,
\]

and the vector part is

\[
\operatorname{Im}(q)=(x,y,z).
\]

Python arrays must use the following order:

```python
q = np.array([w, x, y, z], dtype=float)
```

Scalar-last representations such as

```python
q = np.array([x, y, z, w])
```

must not be used internally without an explicit conversion.

---

## 3. Algebra convention

We use Hamilton quaternion multiplication:

\[
\mathbf i^2=\mathbf j^2=\mathbf k^2=-1,
\]

\[
\mathbf i\mathbf j=\mathbf k,\qquad
\mathbf j\mathbf k=\mathbf i,\qquad
\mathbf k\mathbf i=\mathbf j.
\]

The reversed products have opposite signs.

For

\[
q_1=(a,\mathbf u),\qquad
q_2=(b,\mathbf v),
\]

the product is

\[
q_1q_2=
\left(
ab-\mathbf u\cdot\mathbf v,\,
a\mathbf v+b\mathbf u+\mathbf u\times\mathbf v
\right).
\]

The cross product uses the standard right-handed orientation of \(\mathbb R^3\).

Quaternion multiplication is associative and generally noncommutative.

---

## 4. Conjugation and inversion

Quaternionic conjugation is defined by

\[
q^*=(w,-x,-y,-z).
\]

The norm is

\[
|q|=\sqrt{w^2+x^2+y^2+z^2}.
\]

For every nonzero quaternion,

\[
q^{-1}=\frac{q^*}{|q|^2}.
\]

For a normalized quaternion,

\[
q^{-1}=q^*.
\]

The following operations must remain distinct:

| Operation | Definition |
|---|---|
| Conjugation | \(q^*=(w,-x,-y,-z)\) |
| Negation | \(-q=(-w,-x,-y,-z)\) |
| Inverse | \(q^{-1}=q^*/|q|^2\) |
| Unit inverse | \(q^{-1}=q^*\), if \(|q|=1\) |
| Antipodal inverse lift | \(-q^*\) |

In particular,

\[
q q^*=1,
\]

whereas for unit \(q\),

\[
q(-q^*)=-1.
\]

The antipodal inverse lift represents the same inverse rotation in SO(3), but it is not the inverse element in SU(2).

---

## 5. Rotation convention

We use active rotations acting on column vectors.

A vector

\[
\mathbf r=(r_x,r_y,r_z)
\]

is embedded into the quaternion algebra as

\[
p=(0,r_x,r_y,r_z).
\]

For a unit quaternion \(q\), the rotated vector is obtained from

\[
\boxed{p'=qpq^{-1}.}
\]

The resulting pure quaternion has the form

\[
p'=(0,\mathbf r').
\]

The vector \(\mathbf r'\) is the actively rotated vector.

We use a right-handed Cartesian coordinate system.

---

## 6. Axis-angle convention

For a unit rotation axis

\[
\mathbf n=(n_x,n_y,n_z)
\]

and a right-handed positive angle \(\theta\), the corresponding quaternion is

\[
q=
\left(
\cos\frac{\theta}{2},
n_x\sin\frac{\theta}{2},
n_y\sin\frac{\theta}{2},
n_z\sin\frac{\theta}{2}
\right).
\]

The axis must satisfy

\[
\|\mathbf n\|=1.
\]

The angle is measured in radians.

A zero-angle rotation is represented by

\[
q=(1,0,0,0).
\]

A rotation of \(2\pi\) is represented in SU(2) by

\[
q=(-1,0,0,0).
\]

Although both represent the identity rotation in SO(3), they are distinct elements of SU(2).

---

## 7. Composition order

Quaternion products act from right to left on vectors.

For

\[
q=q_1q_2,
\]

the action is

\[
p'=q_1(q_2pq_2^{-1})q_1^{-1}.
\]

Therefore, \(q_2\) acts first and \(q_1\) acts second.

The rotation matrices obey

\[
R(q_1q_2)=R(q_1)R(q_2).
\]

All implementations must test this identity explicitly.

---

## 8. Rotation matrix convention

The rotation matrix acts on column vectors:

\[
\mathbf r'=R(q)\mathbf r.
\]

For a normalized quaternion \(q=(w,x,y,z)\),

\[
R(q)=
\begin{pmatrix}
1-2(y^2+z^2)&2(xy-wz)&2(xz+wy)\\
2(xy+wz)&1-2(x^2+z^2)&2(yz-wx)\\
2(xz-wy)&2(yz+wx)&1-2(x^2+y^2)
\end{pmatrix}.
\]

The expected properties are

\[
R^TR=I,
\]

\[
\det R=+1.
\]

For every unit quaternion,

\[
R(q)=R(-q).
\]

Implementations must reject or explicitly normalize non-unit inputs rather than silently assuming unit norm.

---

## 9. SU(2) matrix convention

We use the Pauli matrices

\[
\sigma_x=
\begin{pmatrix}
0&1\\
1&0
\end{pmatrix},
\quad
\sigma_y=
\begin{pmatrix}
0&-i\\
i&0
\end{pmatrix},
\quad
\sigma_z=
\begin{pmatrix}
1&0\\
0&-1
\end{pmatrix}.
\]

The canonical quaternion-to-SU(2) map is

\[
U(q)=wI-i(x\sigma_x+y\sigma_y+z\sigma_z).
\]

Explicitly,

\[
U(q)=
\begin{pmatrix}
w-iz&-y-ix\\
y-ix&w+iz
\end{pmatrix}.
\]

For unit \(q\),

\[
U^\dagger U=I,\qquad \det U=1.
\]

The map must preserve multiplication:

\[
U(q_1q_2)=U(q_1)U(q_2).
\]

The antipodal relation is

\[
U(-q)=-U(q).
\]

The two matrices are distinct in SU(2), although they induce the same SO(3) rotation.

---

## 10. Link orientation convention

For an oriented edge from vertex \(g\) to vertex \(h\), define

\[
U_{gh}\in SU(2).
\]

The reverse edge must satisfy

\[
\boxed{U_{hg}=U_{gh}^{-1}=U_{gh}^{\dagger}.}
\]

For a unit quaternion representation,

\[
q_{hg}=q_{gh}^*.
\]

The antipodal representative \(-q_{gh}^*\) must not be substituted for the inverse unless the calculation explicitly tracks the corresponding central sign.

---

## 11. Numerical tolerances

Tests should use scale-aware tolerances.

For standard double-precision computations, an initial target is:

\[
\varepsilon_{\mathrm{unit}}\leq 10^{-12},
\]

\[
\varepsilon_{\mathrm{orth}}\leq 10^{-12},
\]

\[
\varepsilon_{\mathrm{det}}\leq 10^{-12}.
\]

These are proposed software test thresholds, not universal mathematical bounds.

Tests involving long products or ill-conditioned operations must report their actual residuals and may require separately justified tolerances.

---

## 12. Mandatory convention tests

Every implementation must verify:

1. Scalar-first component ordering.
2. Hamilton multiplication.
3. Associativity within numerical precision.
4. Quaternion conjugation.
5. Quaternion inverse.
6. Unit norm.
7. Active rotation of column vectors.
8. Right-handed axis-angle convention.
9. Right-to-left composition order.
10. Rotation matrix orthogonality and determinant.
11. SU(2) unitarity and determinant.
12. The identity \(R(q)=R(-q)\).
13. The distinction between \(q^{-1}\) and \(-q^{-1}\).
14. Reverse-edge inverse consistency.

---

## 13. Scientific status

The conventions in this document define the software interface and mathematical representation used by this repository.

They do not establish any physical interpretation of the broader RGS–PSU or QLG framework.

Alternative conventions are welcome if they are mathematically equivalent, explicitly documented and independently tested.

**Canonical convention:** scalar-first Hamilton quaternions, active right-handed rotations, column vectors, and right-to-left composition.
