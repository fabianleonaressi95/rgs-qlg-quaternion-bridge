# Explicit binary icosahedral group \(2I\)

## Construction

This repository uses scalar-first unit quaternions \(q=(w,x,y,z)\), the
Hamilton product, and active rotations acting on column vectors.

The 120 elements are the vertices of the 600-cell in \(\mathbb R^4\):

1. 8 coordinate units: permutations of \((\pm1,0,0,0)\);
2. 16 half-coordinate points \((\pm\tfrac12,\pm\tfrac12,\pm\tfrac12,\pm\tfrac12)\);
3. 96 points obtained by even permutations of
   \((0,\tfrac12,\tfrac\varphi2,\tfrac{1}{2\varphi})\), with independent
   signs on the three nonzero coordinates, where
   \(\varphi=(1+\sqrt5)/2\).

All have unit norm because
\[
\frac{1+\varphi^2+\varphi^{-2}}4=1.
\]

The implementation is in [src/binary_icosahedral.py](../src/binary_icosahedral.py).
The complete deterministic table is generated from the exact defining pattern;
no numerical search or inferred list is used.

## Exhaustive finite-group audit

The audit checks all \(120^2=14,400\) Hamilton products against the generated
set, verifies conjugate inverses, unit norms, and checks the induced rotation
matrices for orthogonality and determinant +1. It also verifies that the
quaternion-to-rotation map has exactly 60 distinct images: \(q\) and \(-q\)
represent the same element of \(SO(3)\).

This establishes computationally that the listed finite set is closed under
the implemented multiplication and maps two-to-one onto a 60-element rotation
group. The mathematical identification as the binary icosahedral group is
supported by the standard 600-cell construction; the code is an exhaustive
floating-point implementation audit, not a formal proof assistant certificate.

## Reproduction

```bash
python -c "from src.binary_icosahedral import audit_binary_icosahedral; print(audit_binary_icosahedral())"
pytest -q
```

## Boundary of the result

This supplies the explicit \(2I\) quaternion set and verifies its finite-group
and double-cover properties. It does **not** by itself establish that a
particular RGS vierbein, Q6 link field, torsion tensor, or QLG dynamics realizes
this group as its holonomy. That requires a separate representation/connection
audit against the recovered physical matrices.
