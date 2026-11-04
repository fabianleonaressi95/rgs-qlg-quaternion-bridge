# Q6-D5: Maxwell-normalized kinetic term and Liouville measure

**Status:** the finite Q6 spectral sector is inherited from the verified cochain
construction. The electric coefficient below is derived from an explicit
Maxwell normalization and Legendre transform; it is not claimed to be generated
by the spatial-only finite Q6 Wilson action.

## 1. Legendre transform

The compact spatial Wilson action has quadratic term

    S_B = (beta/2) sum_p F_p^2 + O(F^4).

It fixes the magnetic coefficient but contains no time links and therefore
does not determine an electric kinetic term by itself.

Adopt the continuous-time Maxwell normalization

    L = (E^2 - B^2)/(2 g^2),       beta = 1/g^2.

The canonical electric momentum is Pi = E/g^2 = beta E. The Legendre
transform gives

    H = (g^2/2) Pi^2 + B^2/(2g^2)
      = Pi^2/(2 beta) + beta B^2/2.

Therefore the electric Hessian is 1/beta, not an independently chosen unit
coefficient. This is conditional on the declared Maxwell normalization.

## 2. Restriction to the Q6 five-mode eigenspace

The finite cochain audit gives the Coulomb-gauge tangent ker(D0^T), dimension
129, and the isolated top coexact eigenspace

    D1^T D1 v = 12 v,       dim(V5) = 5.

Choose an orthonormal basis V in R^(192x5), with V^T V=I and D0^T V=0.
Write theta=Vq and the transverse canonical momentum Pi=Vp. The reduced
Hamiltonian is

    H5(q,p) = p^T p/(2 beta) + 12 beta q^T q/2.

Its frequency is

    omega5 = sqrt((12 beta)(1/beta)) = sqrt(12),

independent of beta. This differs from the earlier unit-kinetic convention,
where omega=sqrt(12 beta). The difference is a normalization choice, not a
numerical discrepancy.

## 3. Symplectic and Kähler structures

Set a=12 beta, b=1/beta, omega=sqrt(a b). With

    Omega = [[0,I],[-I,0]]
    J = [[0,-(b/omega)I],[(a/omega)I,0]]
    g = Omega J = diag((a/omega)I,(b/omega)I),

the implementation checks J^2=-I, J^T Omega J=Omega and g positive.
The canonical complex coordinate is

    z = (sqrt(a/omega) q + i sqrt(b/omega) p)/sqrt(2),
    H5 = omega z^dagger z.

For the type-IV D5 origin Bergman metric convention g_i_bar_j(0)=10 delta,
the normalized tangent chart is z_D=z/sqrt(10).

## 4. Liouville Jacobian

The real determinant of (q,p) -> (Re z_D, Im z_D) is

    det = 1/(32 * 10^5).

It is independent of beta and of the Hodge eigenvalue because
(a/omega)(b/omega)=1. Thus

    d^5q d^5p = 32 * 10^5 d^10(Re z_D, Im z_D).

This is only a constant Jacobian for the linear tangent chart. It does not
show that the compact 129-dimensional quotient Haar measure pushes forward
to the global Bergman measure on noncompact D5. No global nonlinear map has
been constructed.

## 5. Evidence ledger

**Derived, conditional on the declared Maxwell normalization**
- electric coefficient 1/(2 beta);
- reduced five-mode Hamiltonian;
- beta-independent frequency sqrt(12);
- compatible complex/Kähler structure;
- constant linear Liouville Jacobian.

**Numerically verified**
- inherited Q6 five-mode basis is orthonormal and Coulomb transverse;
- symplectic, complex and Kähler identities;
- energy identity in complex coordinates;
- Jacobian for multiple positive beta values.

**Still open**
- deriving Maxwell normalization and beta from the RGS-PSU action;
- adding temporal links/plaquettes and deriving electric energy from finite Q6;
- global Q6-to-D5 embedding and global measure push-forward;
- continuum limit and matching to the Thomson coupling / fine-structure constant.

## Reproduction

    PYTHONPATH=. python -m pytest tests/test_rgs_q6_d5_maxwell_hamiltonian.py -q
    PYTHONPATH=. python -c "from src.rgs_q6_d5_maxwell_hamiltonian import build_q6_d5_maxwell_hamiltonian; print(build_q6_d5_maxwell_hamiltonian()['all_checks_passed'])"
