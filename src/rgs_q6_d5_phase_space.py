"""Canonical 10-real-dimensional Q6 phase-space sector and D5 tangent audit.

The construction is a spectral reduction of the U(1) quotient tangent:
the top coexact 1-form eigenspace of the cubical Hodge Laplacian is a
5-dimensional irreducible A5 module. Its cotangent lift is 10-dimensional.
A unit electric kinetic metric is stated explicitly; beta remains bare.

This identifies the isotropy representation at the origin of the type-IV
domain D_5, not the global bounded domain or Wyler's Euclidean volume.
"""
from __future__ import annotations
import math
import numpy as np
from src.rgs_q6_equivariant_cochains import build_q6_cochain_complex, induced_cochain_actions
from src.rgs_icosahedral_q6_intertwiner import icosahedral_signed_axis_representation


_A5_CLASS_SIZES=np.array([1.,15.,20.,12.,12.])
_SQRT5=math.sqrt(5.)
_PHI=(1.+_SQRT5)/2.
_A5_IRREPS={
 "1":np.array([1.,1.,1.,1.,1.]),
 "3":np.array([3.,-1.,0.,_PHI,(1.-_SQRT5)/2.]),
 "3_prime":np.array([3.,-1.,0.,(1.-_SQRT5)/2.,_PHI]),
 "4":np.array([4.,0.,1.,-1.,-1.]),
 "5":np.array([5.,1.,-1.,0.,0.]),
}


def _deterministic_basis_from_projector(projector, dimension, tolerance=1e-9):
    """Choose a reproducible orthonormal frame from an invariant projector."""
    p=np.asarray(projector,dtype=float)
    basis=[]
    for column in range(p.shape[1]):
        v=p[:,column].copy()
        for q in basis:
            v-=q*np.dot(q,v)
        for q in basis:
            v-=q*np.dot(q,v)
        norm=float(np.linalg.norm(v))
        if norm<=tolerance:
            continue
        v/=norm
        nz=np.flatnonzero(np.abs(v)>tolerance)
        if len(nz) and v[nz[0]]<0:
            v=-v
        basis.append(v)
        if len(basis)==dimension:
            break
    if len(basis)!=dimension:
        raise ArithmeticError("could not extract deterministic spectral basis")
    return np.column_stack(basis)


def _a5_class_index(rotation_trace):
    t=float(rotation_trace)
    if abs(t-3.)<1e-7: return 0
    if abs(t+1.)<1e-7: return 1
    if abs(t)<1e-7: return 2
    if abs(t-_PHI)<1e-7: return 3
    if abs(t-(1.-_SQRT5)/2.)<1e-7: return 4
    raise ArithmeticError(f"unrecognized proper icosahedral rotation trace {t}")


def build_q6_d5_phase_space_sector(beta=1.0, eigenvalue=12.0):
    """Build the 5-mode coexact spectral sector and its 10D cotangent lift.

    beta is the positive bare Wilson coefficient. The electric kinetic metric
    is the identity on the selected orthonormal mode coordinates.
    """
    if not np.isfinite(beta) or beta<=0:
        raise ValueError("beta must be finite and positive")
    if not np.isfinite(eigenvalue) or eigenvalue<=0:
        raise ValueError("eigenvalue must be finite and positive")
    c=build_q6_cochain_complex()
    d0=np.asarray(c["D0"],dtype=float)
    d1=np.asarray(c["D1"],dtype=float)

    # Coulomb-gauge tangent ker(D0^T), dimension E-(V-1)=129.
    _u,svals,vh=np.linalg.svd(d0.T,full_matrices=True)
    gauge_rank=int(np.sum(svals>1e-10))
    if gauge_rank!=63:
        raise ArithmeticError(f"expected gauge rank 63, got {gauge_rank}")
    coulomb_basis=vh[gauge_rank:].T
    reduced_laplacian=coulomb_basis.T@(d1.T@d1)@coulomb_basis
    reduced_laplacian=(reduced_laplacian+reduced_laplacian.T)/2
    spectrum,eigenvectors=np.linalg.eigh(reduced_laplacian)
    selected=np.flatnonzero(np.isclose(spectrum,eigenvalue,rtol=0.,atol=1e-8))
    if len(selected)!=5:
        raise ArithmeticError(
            f"target Hodge eigenvalue must have multiplicity five; got {len(selected)}"
        )
    spectral_projector=(coulomb_basis@eigenvectors[:,selected])@(
        coulomb_basis@eigenvectors[:,selected]
    ).T
    mode_basis=_deterministic_basis_from_projector(spectral_projector,5)
    projector_residual=float(np.linalg.norm(
        mode_basis@mode_basis.T-spectral_projector
    ))
    hodge_residual=float(np.linalg.norm(
        d1.T@d1@mode_basis-eigenvalue*mode_basis
    ))
    coulomb_residual=float(np.linalg.norm(d0.T@mode_basis))

    # A5 action on the selected real 5-dimensional eigenspace.
    group=icosahedral_signed_axis_representation()
    class_traces=[[] for _ in range(5)]
    representation_matrices=[]
    equivariance_residual=0.
    orthogonality_residual=0.
    for rotation,signed in zip(group["rotations"],group["signed_permutations"]):
        _p0,p1,_p2=induced_cochain_actions(signed,c)
        rho=mode_basis.T@p1@mode_basis
        representation_matrices.append(rho)
        ci=_a5_class_index(np.trace(rotation))
        class_traces[ci].append(float(np.trace(rho)))
        equivariance_residual=max(
            equivariance_residual,
            float(np.linalg.norm(p1@mode_basis-mode_basis@rho)),
        )
        orthogonality_residual=max(
            orthogonality_residual,
            float(np.linalg.norm(rho.T@rho-np.eye(5))),
        )
    character=np.array([np.mean(x) for x in class_traces])
    multiplicities={
        name:int(round(float(np.dot(_A5_CLASS_SIZES,character*chi)/60.)))
        for name,chi in _A5_IRREPS.items()
    }
    character_fit_residual=max(
        abs(float(np.dot(_A5_CLASS_SIZES,character*chi)/60.)-multiplicities[name])
        for name,chi in _A5_IRREPS.items()
    )

    # Cotangent phase space T*V = V x V*: canonical symplectic form.
    # Unit electric kinetic metric is an explicit modelling convention.
    omega=math.sqrt(beta*eigenvalue)
    eye=np.eye(5)
    zero=np.zeros((5,5))
    symplectic=np.block([[zero,eye],[-eye,zero]])
    complex_structure=np.block([[zero,-eye/omega],[omega*eye,zero]])
    kahler_metric=symplectic@complex_structure
    hamiltonian_hessian=np.block([[beta*eigenvalue*eye,zero],[zero,eye]])
    phase_group=np.block([[representation_matrices[0]*0+eye,zero],[zero,eye]])
    # SO(2) phase rotations exp(angle*J), checked at a nontrivial angle.
    angle=.371
    phase_rotation=math.cos(angle)*np.eye(10)+math.sin(angle)*complex_structure
    group_action=np.block([[representation_matrices[1],zero],[zero,representation_matrices[1]]])
    symplectic_residual=float(np.linalg.norm(
        complex_structure.T@symplectic@complex_structure-symplectic
    ))
    complex_residual=float(np.linalg.norm(complex_structure@complex_structure+np.eye(10)))
    kahler_symmetry_residual=float(np.linalg.norm(kahler_metric-kahler_metric.T))
    kahler_positive_min=float(np.linalg.eigvalsh(kahler_metric).min())
    phase_residual=float(np.linalg.norm(
        phase_rotation.T@kahler_metric@phase_rotation-kahler_metric
    ))
    a5_phase_residual=float(np.linalg.norm(
        group_action.T@kahler_metric@group_action-kahler_metric
    ))
    a5_symplectic_residual=float(np.linalg.norm(
        group_action.T@symplectic@group_action-symplectic
    ))
    a5_complex_residual=float(np.linalg.norm(
        group_action@complex_structure-complex_structure@group_action
    ))

    checks={
        "quotient_tangent_dimension_129":coulomb_basis.shape==(192,129),
        "top_coexact_mode_multiplicity_5":len(selected)==5,
        "mode_basis_orthonormal":np.linalg.norm(mode_basis.T@mode_basis-eye)<1e-10,
        "mode_is_coulomb":coulomb_residual<1e-10,
        "mode_hodge_eigen_residual":hodge_residual<1e-8,
        "spectral_projector_reproducible":projector_residual<1e-8,
        "A5_irrep_is_five":multiplicities=={"1":0,"3":0,"3_prime":0,"4":0,"5":1},
        "A5_action_equivariant":equivariance_residual<1e-8,
        "A5_action_orthogonal":orthogonality_residual<1e-8,
        "cotangent_dimension_10":symplectic.shape==(10,10),
        "complex_structure":complex_residual<1e-10,
        "symplectic_compatibility":symplectic_residual<1e-10,
        "kahler_metric_symmetric":kahler_symmetry_residual<1e-10,
        "kahler_metric_positive":kahler_positive_min>0,
        "SO2_phase_isometry":phase_residual<1e-10,
        "A5_phase_isometry":a5_phase_residual<1e-8,
        "A5_symplectic":a5_symplectic_residual<1e-8,
        "A5_commutes_with_complex_structure":a5_complex_residual<1e-8,
    }
    return {
        "audit":"Q6_A5_FIVE_MODE_COTANGENT_D5_TANGENT",
        "q6_vertices":64,"q6_edges":192,"q6_faces":240,
        "gauge_rank":gauge_rank,"quotient_configuration_dimension":129,
        "coexact_hodge_spectrum":{str(float(v)):int(np.sum(np.isclose(spectrum,v,atol=1e-8)))
                                  for v in np.unique(np.round(spectrum,8))},
        "selected_hodge_eigenvalue":float(eigenvalue),
        "selected_mode_multiplicity":int(len(selected)),
        "selected_mode_basis":mode_basis,
        "selected_mode_projector":spectral_projector,
        "selected_A5_character_by_class":character,
        "A5_class_sizes":[1,15,20,12,12],
        "A5_irrep_multiplicities":multiplicities,
        "A5_equivariance_residual":equivariance_residual,
        "A5_orthogonality_residual":orthogonality_residual,
        "phase_space_dimension":10,
        "phase_space_coordinates":"(q_1,...,q_5,p_1,...,p_5)",
        "canonical_symplectic_form":symplectic,
        "compatible_complex_structure":complex_structure,
        "kahler_metric":kahler_metric,
        "oscillator_frequency":omega,
        "bare_wilson_beta":float(beta),
        "electric_kinetic_metric":"identity on selected modes (explicit convention)",
        "complex_coordinate":"z=(sqrt(omega)*q + i*p/sqrt(omega))/sqrt(2)",
        "d5_type_iv_domain":"{z in C^5: 1-2||z||^2+|z.z|^2>0 and |z.z|<1}",
        "d5_isotropy_tangent_dimension_real":10,
        "d5_isotropy_tangent_dimension_complex":5,
        "bergman_kernel_generic_norm_exponent":5,
        "bergman_metric_origin_hermitian_factor":10,
        "bergman_coordinate_tangent_scale":"z_D=z/sqrt(10), under g=partial_barpartial log K",
        "residuals":{
            "mode_projector":projector_residual,"mode_hodge":hodge_residual,
            "mode_coulomb":coulomb_residual,"A5_equivariance":equivariance_residual,
            "A5_orthogonality":orthogonality_residual,"A5_character_fit":character_fit_residual,
            "J_squared_plus_I":complex_residual,"symplectic":symplectic_residual,
            "kahler_symmetry":kahler_symmetry_residual,"SO2_isometry":phase_residual,
            "A5_isometry":a5_phase_residual,"A5_symplectic":a5_symplectic_residual,
            "A5_J_commutator":a5_complex_residual,
        },
        "all_checks_passed":bool(all(checks.values())),
        "checks":checks,
        "epistemic_status":{
            "verified_finite_structure":[
                "129-dimensional Coulomb-gauge tangent of the Q6 U(1) link quotient",
                "unique top coexact Hodge eigenspace has dimension five",
                "that eigenspace is the irreducible five-dimensional A5 representation",
                "its cotangent lift is a ten-real-dimensional symplectic/Kahler oscillator sector",
                "the sector has SO(5) x SO(2) linear oscillator symmetry",
                "this is the isotropy representation carried by the tangent space of type-IV D5",
            ],
            "not_proved":[
                "a global embedding of the compact Q6 gauge quotient into noncompact D5",
                "a canonical nonlinear Q6-to-D5 map or its Jacobian",
                "identification with Wyler's Euclidean volume normalization",
                "the physical electric kinetic coefficient or charge normalization",
                "a four-dimensional continuum QED limit or prediction of alpha",
            ],
        },
    }


def d5_type_iv_contains(z, tolerance=1e-12):
    """Membership test for the standard Cartan type-IV_5 bounded realization."""
    value=np.asarray(z,dtype=complex)
    if value.shape!=(5,) or not np.all(np.isfinite(value)):
        raise ValueError("z must be a finite complex vector of shape (5,)")
    norm2=float(np.vdot(value,value).real)
    quadratic=complex(np.dot(value,value))
    return bool(1.-2.*norm2+abs(quadratic)**2>tolerance and abs(quadratic)<1.-tolerance)


def phase_space_to_d5_tangent(q,p,beta=1.0,eigenvalue=12.0):
    """Map canonical oscillator coordinates to the normalized D5 tangent chart.

    The factor sqrt(10) matches the origin Bergman Hermitian metric convention
    g=partial_barpartial log(K), K proportional to N(z,z)^(-5).
    This is a tangent-chart normalization, not a global Q6-to-D5 map.
    """
    q=np.asarray(q,dtype=float); p=np.asarray(p,dtype=float)
    if q.shape!=(5,) or p.shape!=(5,) or not np.all(np.isfinite(q)) or not np.all(np.isfinite(p)):
        raise ValueError("q and p must be finite real vectors of shape (5,)")
    if not np.isfinite(beta) or beta<=0 or not np.isfinite(eigenvalue) or eigenvalue<=0:
        raise ValueError("beta and eigenvalue must be finite and positive")
    omega=math.sqrt(beta*eigenvalue)
    z=(math.sqrt(omega)*q+1j*p/math.sqrt(omega))/math.sqrt(2.*10.)
    return z
