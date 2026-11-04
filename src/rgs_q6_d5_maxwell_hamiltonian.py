from __future__ import annotations
import math
import numpy as np
from src.rgs_q6_d5_phase_space import build_q6_d5_phase_space_sector, d5_type_iv_contains
from src.rgs_q6_cochain_complex import build_q6_cochain_complex


def maxwell_mode_parameters(beta=1.0, eigenvalue=12.0):
    beta, eigenvalue = float(beta), float(eigenvalue)
    if not math.isfinite(beta) or beta <= 0 or not math.isfinite(eigenvalue) or eigenvalue <= 0:
        raise ValueError('beta and eigenvalue must be finite and positive')
    electric = 1.0 / beta
    magnetic = beta * eigenvalue
    omega = math.sqrt(electric * magnetic)
    return {'beta': beta, 'g_squared': 1.0 / beta, 'electric_hessian': electric,
            'magnetic_hessian': magnetic, 'oscillator_frequency': omega,
            'q_scale': math.sqrt(magnetic / omega / 2.0),
            'p_scale': math.sqrt(electric / omega / 2.0),
            'kinetic_coefficient': 1.0 / (2.0 * beta),
            'magnetic_coefficient': beta * eigenvalue / 2.0,
            'normalization': 'beta=1/g^2; canonical Pi=E/g^2'}


def phase_space_to_d5_maxwell_tangent(q, p, beta=1.0, eigenvalue=12.0):
    q, p = np.asarray(q, float), np.asarray(p, float)
    if q.shape != (5,) or p.shape != (5,) or not np.all(np.isfinite(q)) or not np.all(np.isfinite(p)):
        raise ValueError('q and p must be finite vectors of shape (5,)')
    par = maxwell_mode_parameters(beta, eigenvalue)
    return (par['q_scale'] * q + 1j * par['p_scale'] * p) / math.sqrt(10.0)


def build_q6_d5_maxwell_hamiltonian(beta=1.0, eigenvalue=12.0):
    par = maxwell_mode_parameters(beta, eigenvalue)
    finite = build_q6_d5_phase_space_sector(beta=beta, eigenvalue=eigenvalue)
    basis = np.asarray(finite['selected_mode_basis'], float)
    d0 = np.asarray(build_q6_cochain_complex()['D0'], float)
    eye, zero = np.eye(5), np.zeros((5, 5))
    omega = par['oscillator_frequency']
    a, b = par['magnetic_hessian'], par['electric_hessian']
    symplectic = np.block([[zero, eye], [-eye, zero]])
    complex_structure = np.block([[zero, -(b / omega) * eye], [(a / omega) * eye, zero]])
    kahler = symplectic @ complex_structure
    q = np.array([.07, -.03, .02, .01, -.04]); p = np.array([-.02, .05, .01, -.06, .03])
    z = phase_space_to_d5_maxwell_tangent(q, p, beta, eigenvalue)
    h = .5 * (a * np.dot(q, q) + b * np.dot(p, p))
    hz = omega * 10.0 * np.vdot(z, z).real
    jac = float(np.linalg.det(np.diag([par['q_scale']/math.sqrt(10.0)]*5 + [par['p_scale']/math.sqrt(10.0)]*5)))
    expected_jac = 1.0 / (32.0 * 10.0**5)
    checks = {
        'basis_orthonormal': np.linalg.norm(basis.T @ basis-eye) < 1e-10,
        'gauss_coulomb': np.linalg.norm(d0.T @ basis) < 1e-9,
        'electric_inverse_beta': math.isclose(b, 1.0/beta, rel_tol=1e-14),
        'magnetic_wilson': math.isclose(a, beta*eigenvalue, rel_tol=1e-14),
        'frequency_beta_independent': math.isclose(omega, math.sqrt(eigenvalue), rel_tol=1e-14),
        'complex_structure': np.linalg.norm(complex_structure@complex_structure+np.eye(10)) < 1e-10,
        'symplectic': np.linalg.norm(complex_structure.T@symplectic@complex_structure-symplectic) < 1e-10,
        'kahler_symmetric': np.linalg.norm(kahler-kahler.T) < 1e-10,
        'kahler_positive': np.linalg.eigvalsh(kahler).min() > 0,
        'liouville_jacobian': math.isclose(jac, expected_jac, rel_tol=1e-12, abs_tol=1e-16),
        'energy_identity': math.isclose(h, hz, rel_tol=1e-12, abs_tol=1e-14),
        'D5_local_chart': d5_type_iv_contains(z),
        'finite_Q6_audit': bool(finite['all_checks_passed'])}
    return {'audit': 'Q6_D5_MAXWELL_NORMALIZED_HAMILTONIAN', 'parameters': par,
            'mode_basis': basis, 'hamiltonian_hessian': np.block([[a*eye,zero],[zero,b*eye]]),
            'symplectic_form': symplectic, 'complex_structure': complex_structure,
            'kahler_metric': kahler, 'liouville_jacobian': jac,
            'expected_liouville_jacobian': expected_jac, 'sample_hamiltonian': h,
            'sample_complex_energy': hz, 'checks': checks,
            'all_checks_passed': all(checks.values()),
            'scope': 'Maxwell Legendre transform conditional on beta=1/g^2; temporal plaquette dynamics is not in the spatial-only Q6 Wilson module.'}
