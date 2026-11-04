#!/usr/bin/env python3
"""
RGS-PSU Closure Script: Mathematical Structure and Geometric Background Solver

STATUS: CORRECTED - DOCUMENTING INCONSISTENCIES

This script implements the mathematical framework for RGS-PSU background analysis
and documents the four critical inconsistencies identified in the implementation.

CORRECTIONS APPLIED:
A. Teleparallel identity: Original formulas for T and ∇_μT^μ were algebraically
   incompatible with R_LC. The identity now correctly yields R_LC + T - 2∇_μT^μ = -2R''/R,
   which only vanishes for special cases, not general R(u).

B. Tetrad correction: For metric ds² = -N(u)²dt² + du² + R(u)²dΩ², the correct
   diagonal tetrad is e^a_μ = diag(N(u), 1, R(u), R(u)sinθ), not diag(1,1,R,Rsinθ).

C. Electromagnetic field: The density ρ_EM = Q²/(2√(r⁴+ℓ₀⁴)) recovers Coulomb ~Q²/(2r²),
   but the reported field E(r) = Qr²/(r⁴+ℓ₀⁴)^(3/2) ~ Q/r⁴ does NOT recover Q/r².
   These cannot both come from the same NLED Lagrangian without derivation.

D. Misner-Sharp mass: For R(0)=2.31 and R'(0)=0, the correct m_MS(0) = R(0)/2 = 1.155,
   not 5.124. The previous value came from an inconsistent implementation.

OPEN PROBLEMS:
- Complete teleparallel gravitational equation for N(u) dynamics
- Derive NLED Lagrangian that yields consistent ρ_EM and E(r)
- Verify which of the three R(0) configurations correspond to the same physical background
- Complete BVP solver for coupled R(u), φ(u), N(u) system

THREE BACKGROUNDS (DOCUMENTED AS SEPARATE):
1. R(0)=2.31, φ(0)=0.80 - Percorso A BVP result (FORMULAZIONE.odt)
2. R(0)=0.8292, N(0)=1.4529, φ(0)=0.3450 - Reported in FORMULAZIONE.pdf (missing)
3. R(0)=2, N(0)=0.9821, φ(0)=0.3450 - Reported in siiiiii.pdf (needs verification)

Author: RGS-PSU-PRISM Program
Date: 2026-10-28
Status: CORRECTED - Inconsistencies documented, not solved
"""

import numpy as np
import json
from scipy.integrate import solve_ivp, quad, solve_bvp
from scipy.optimize import minimize, root
from scipy.linalg import norm
import warnings
warnings.filterwarnings('ignore')

# ==============================================================================
# SECTION 1: GEOMETRIC IDENTITIES (EXACT, INDEPENDENT OF ACTION)
# ==============================================================================

class GeometricIdentities:
    """
    Exact geometric identities for static spherically symmetric metric:
    ds^2 = -N(u)^2 dt^2 + du^2 + R(u)^2 dΩ^2
    
    Includes both Einstein tensor components AND teleparallel torsion quantities
    from FORMULAZIONE.odt (Percorso A).
    """
    
    @staticmethod
    def G_tt(N, R, u):
        """G^t_t = 2R''/R + R'^2/R^2 - 1/R^2"""
        Rp = np.gradient(R, u)
        Rpp = np.gradient(Rp, u)
        return 2*Rpp/R + Rp**2/R**2 - 1/R**2
    
    @staticmethod
    def G_uu(N, R, u):
        """G^u_u = 2N'R'/(NR) + R'^2/R^2 - 1/R^2"""
        Rp = np.gradient(R, u)
        Np = np.gradient(N, u)
        return 2*Np*Rp/(N*R) + Rp**2/R**2 - 1/R**2
    
    @staticmethod
    def G_thetatheta(N, R, u):
        """G^θ_θ = N''/N + N'R'/(NR) + R''/R"""
        Rp = np.gradient(R, u)
        Np = np.gradient(N, u)
        Npp = np.gradient(Np, u)
        Rpp = np.gradient(Rp, u)
        return Npp/N + Np*Rp/(N*R) + Rpp/R
    
    @staticmethod
    def torsion_scalar(R, u):
        """
        Teleparallel torsion scalar T(u) from FORMULAZIONE.odt.
        
        For tetrad e^a_μ = diag(1, 1, R(u), R(u)sinθ):
        T(u) = -2R''/R - 2R'²/R² - 4/R²
        """
        Rp = np.gradient(R, u)
        Rpp = np.gradient(Rp, u)
        return -2*Rpp/R - 2*Rp**2/R**2 - 4/R**2
    
    @staticmethod
    def torsion_vector_divergence(R, u):
        """
        Divergence of torsion vector ∇_μT^μ from FORMULAZIONE.odt.
        
        For tetrad e^a_μ = diag(1, 1, R(u), R(u)sinθ):
        ∇_μT^μ = -2R''/R - 2R'²/R² - 1/R²
        """
        Rp = np.gradient(R, u)
        Rpp = np.gradient(Rp, u)
        return -2*Rpp/R - 2*Rp**2/R**2 - 1/R**2
    
    @staticmethod
    def levi_civita_scalar(R, u):
        """
        Levi-Civita curvature scalar R_LC for ultrastatic metric:
        ds² = -dt² + du² + R(u)²dΩ²
        
        CORRECTED FORMULA:
        R_LC = -4R''/R + 2(1 - R'²)/R²
        """
        Rp = np.gradient(R, u)
        Rpp = np.gradient(Rp, u)
        return -4*Rpp/R + 2*(1 - Rp**2)/R**2
    
    @staticmethod
    def verify_teleparallel_identity(R, u):
        """
        Verify the teleparallel identity: R_LC + T - 2∇_μT^μ
        
        CORRECTED: This equals -2R''/R, NOT zero in general.
        Returns the computed value of the identity across the grid.
        """
        R_LC = GeometricIdentities.levi_civita_scalar(R, u)
        T = GeometricIdentities.torsion_scalar(R, u)
        div_T = GeometricIdentities.torsion_vector_divergence(R, u)
        
        identity_value = R_LC + T - 2*div_T
        return {
            "max_absolute": np.max(np.abs(identity_value)),
            "mean_value": np.mean(identity_value),
            "theoretical_value": -2 * np.gradient(np.gradient(R, u), u) / R  # -2R''/R
        }


# ==============================================================================
# SECTION 2: QUANTITIES TO COMPUTE (MODEL-DEPENDENT)
# ==============================================================================

class PhysicalQuantities:
    """
    Physical quantities computed from the background solution.
    These depend on the specific action and matter content.
    """
    
    @staticmethod
    def misner_sharp_mass(R, G=1.0):
        """
        Misner-Sharp mass: m_MS(u) = R/(2G) * (1 - R'^2)
        
        CORRECTED: At throat where R'(0)=0, m_MS(0) = R(0)/2.
        For R(0)=2.31, this gives m_MS(0)=1.155, not 5.124.
        
        Args:
            R: radial profile function
            G: gravitational constant (default 1 in geometric units)
        
        Returns:
            m_MS(u) along the grid
        """
        Rp = np.gradient(R)
        return R / (2*G) * (1 - Rp**2)
    
    @staticmethod
    def electromagnetic_energy_regularized(R, Q, ell0=1.0, G=1.0):
        """
        Electromagnetic energy using RGS-EMS regularized density from siiiiii.pdf:
        ρ_EM(r) = Q² / (2√(r⁴ + ℓ₀⁴))
        
        This is the RGS-EMS regularization (NOT provisional).
        
        Args:
            R: radial profile
            Q: charge parameter
            ell0: regularization scale (core size)
            G: gravitational constant
        
        Returns:
            E_EM integrated over two symmetric branches
        """
        rho_em = Q**2 / (2 * np.sqrt(R**4 + ell0**4))
        integrand = 8 * np.pi * rho_em * R**2
        
        # Integrate over the domain (symmetric)
        u = np.linspace(0, len(R)-1, len(R))
        E_em = np.trapz(integrand, u)
        
        return E_em, rho_em
    
    @staticmethod
    def scalar_potential(phi, V0=1.0, lam=1.0, phi_v=1.0):
        """
        Double-well potential from siiiiii.pdf:
        V(φ) = V₀ + (λ/4)(φ² - φ_v²)²
        
        Derivatives:
        V'(φ) = λφ(φ² - φ_v²)
        V''(φ) = λ(3φ² - φ_v²)
        """
        V = V0 + (lam/4) * (phi**2 - phi_v**2)**2
        V_prime = lam * phi * (phi**2 - phi_v**2)
        V_double_prime = lam * (3*phi**2 - phi_v**2)
        return V, V_prime, V_double_prime
    
    @staticmethod
    def compute_residuals(equations, u):
        """
        Compute L2 residual for each equation.
        
        For each equation E_i(u) = 0, compute:
        Res_i = (integral |E_i(u)|^2 du)^(1/2)
        
        Args:
            equations: dict of equation names to arrays E_i(u)
            u: coordinate grid
        
        Returns:
            dict of residuals
        """
        residuals = {}
        for name, E in equations.items():
            residuals[name] = np.sqrt(np.trapz(np.abs(E)**2, u))
        return residuals


# ==============================================================================
# SECTION 3: ADM MASS EXTRACTION WITH CONVERGENCE TESTS
# ==============================================================================

class ADMExtraction:
    """
    ADM mass extraction from asymptotic behavior of m_MS(u).
    
    Using diagnostic fit: m_MS(R) = M + c1/R + c2/R^2
    The coefficient M is accepted as ADM candidate ONLY after convergence tests.
    """
    
    def __init__(self, convergence_tolerance=1e-6):
        self.tolerance = convergence_tolerance
        self.convergence_report = {}
    
    def asymptotic_fit(self, R, m_MS, R_min=None):
        """
        Fit m_MS(R) = M + c1/R + c2/R^2 in asymptotic region.
        
        Args:
            R: radial coordinate
            m_MS: Misner-Sharp mass
            R_min: minimum R for asymptotic region (auto-determined if None)
        
        Returns:
            dict with fit parameters and convergence metrics
        """
        # Determine asymptotic region where R' ≈ 1
        Rp = np.gradient(R)
        asymptotic_mask = (Rp > 0.95) & (Rp < 1.05)
        
        if R_min is not None:
            asymptotic_mask = asymptotic_mask & (R > R_min)
        
        R_asym = R[asymptotic_mask]
        m_asym = m_MS[asymptotic_mask]
        
        if len(R_asym) < 10:
            return {"status": "insufficient_asymptotic_points", "M_candidate": None}
        
        # Perform fit
        def model(R, M, c1, c2):
            return M + c1/R + c2/R**2
        
        from scipy.optimize import curve_fit
        try:
            popt, pcov = curve_fit(model, R_asym, m_asym, p0=[m_asym[-1], 0, 0])
            M_candidate, c1, c2 = popt
            perr = np.sqrt(np.diag(pcov))
            
            result = {
                "status": "success",
                "M_candidate": M_candidate,
                "M_error": perr[0],
                "c1": c1,
                "c2": c2,
                "n_points": len(R_asym),
                "R_range": [R_asym.min(), R_asym.max()]
            }
            
            return result
        except Exception as e:
            return {"status": f"fit_failed: {str(e)}", "M_candidate": None}
    
    def test_convergence(self, R_profiles, m_MS_profiles):
        """
        Test convergence of ADM mass extraction with respect to:
        - Number of grid points
        - Domain extension
        - Asymptotic behavior R' → 1
        
        Args:
            R_profiles: list of R(u) profiles with different resolutions
            m_MS_profiles: corresponding m_MS(u) profiles
        
        Returns:
            convergence report
        """
        M_candidates = []
        
        for R, m_ms in zip(R_profiles, m_MS_profiles):
            fit_result = self.asymptotic_fit(R, m_ms)
            if fit_result["status"] == "success":
                M_candidates.append(fit_result["M_candidate"])
        
        if len(M_candidates) < 2:
            return {"status": "insufficient_profiles_for_convergence"}
        
        # Check convergence
        M_std = np.std(M_candidates)
        M_mean = np.mean(M_candidates)
        
        convergence = {
            "M_values": M_candidates,
            "M_mean": M_mean,
            "M_std": M_std,
            "relative_std": M_std / abs(M_mean) if M_mean != 0 else np.inf,
            "converged": M_std < self.tolerance * abs(M_mean)
        }
        
        return convergence


# ==============================================================================
# SECTION 4: BACKGROUND SOLVER (SHOOTING METHOD)
# ==============================================================================

class BackgroundSolver:
    """
    Solver for RGS-EMS background equations using actual field equations.
    
    Implements equations from:
    - siiiiii.pdf: RGS-EMS unified action
    - FORMULAZIONE.odt: Percorso A teleparallel scalar-tensor theory
    """
    
    def __init__(self, G=1.0, Q=1.0, ell0=1.0, xi=0.5, V0=1.0, lam=1.0, phi_v=1.0):
        self.G = G
        self.Q = Q
        self.ell0 = ell0
        self.xi = xi  # Scalar-torsion coupling F(φ) = 1 + ξφ²
        self.V0 = V0
        self.lam = lam
        self.phi_v = phi_v
        self.geometry = GeometricIdentities()
        self.quantities = PhysicalQuantities()
        self.adm = ADMExtraction()
    
    def scalar_field_equation(self, R, phi, u):
        """
        Scalar field equation from FORMULAZIONE.odt (Percorso A):
        □φ - V'(φ) - ½F'(φ)T + ∇_μ(F'(φ)T^μ) = 0
        
        For F(φ) = 1 + ξφ²:
        F'(φ) = 2ξφ
        Term becomes: -ξφT + ∇_μ(2ξφT^μ)
        
        In static spherically symmetric coordinates:
        (1/R²)∂_u(R²N∂_uφ) - N V'(φ) - ξφT + 2ξ∇_μ(φT^μ) = 0
        """
        Rp = np.gradient(R, u)
        Rpp = np.gradient(Rp, u)
        phip = np.gradient(phi, u)
        phipp = np.gradient(phip, u)
        
        # Compute torsion quantities
        T = self.geometry.torsion_scalar(R, u)
        div_T = self.geometry.torsion_vector_divergence(R, u)
        
        # F'(φ) for F(φ) = 1 + ξφ²
        F_prime = 2 * self.xi * phi
        
        # Potential derivative
        _, V_prime, _ = self.quantities.scalar_potential(phi, self.V0, self.lam, self.phi_v)
        
        # Simplified: assume N≈1 for placeholder (should be computed from full system)
        N = 1.0
        
        # D'Alembertian term: (1/R²)∂_u(R²∂_uφ)
        dalembert_phi = (1/R**2) * np.gradient(R**2 * phip, u)
        
        # Full scalar equation
        scalar_eq = dalembert_phi - N * V_prime - 0.5 * F_prime * T + self.xi * phi * div_T
        
        return scalar_eq
    
    def field_equations_percorso_A(self, u, y):
        """
        Full system of ODEs for Percorso A from FORMULAZIONE.odt.
        
        State vector y = [R, R', φ, φ', N]
        
        Equations:
        1. R' = y[1]
        2. R'' from teleparallel Hamiltonian constraint
        3. φ' = y[3]
        4. φ'' from scalar equation
        5. N' from Hamiltonian constraint
        
        Action: S = ∫ d⁴x e[(1+ξφ²)T/2 - T^μ∂_μ(1+ξφ²) - ½(∂φ)² - V(φ)]
        
        Gravitational equation (tt component):
        Using G^t_t = 8πG T^t_t with teleparallel torsion corrections
        For metric ds² = -N²dt² + du² + R²dΩ²:
        G^t_t = 2R''/R + R'²/R² - 1/R²
        
        With scalar source: T^t_t = -ρ_φ - ρ_EM
        ρ_φ = (1/2)N²φ'² + V(φ) + torsion coupling terms
        """
        R = y[0]
        Rp = y[1]
        phi = y[2]
        phip = y[3]
        N = y[4]
        
        # Torsion quantities (depend on R and derivatives)
        # Rpp needed for T and div_T - this creates coupling
        # For now, use current approximation then iterate
        Rpp = 0.0  # Will be determined from constraint
        
        # Re-compute torsion with Rpp
        T = -2*Rpp/R - 2*Rp**2/R**2 - 4/R**2
        div_T = -2*Rpp/R - 2*Rp**2/R**2 - 1/R**2
        
        # Scalar potential and derivatives
        V, V_prime, V_double_prime = self.quantities.scalar_potential(phi, self.V0, self.lam, self.phi_v)
        
        # Scalar density (kinetic + potential)
        rho_phi = 0.5 * N**2 * phip**2 + V
        
        # EM density (regularized)
        rho_em = self.Q**2 / (2 * np.sqrt(R**4 + self.ell0**4))
        
        # Total source density
        rho_total = rho_phi + rho_em
        
        # Teleparallel Hamiltonian constraint: G^t_t = 8πG T^t_t
        # 2R''/R + R'²/R² - 1/R² = 8πG (-rho_total)
        # Solve for R'':
        Rpp = (R/2) * (8*np.pi*self.G * rho_total - Rp**2/R**2 + 1/R**2)
        
        # Now re-compute torsion with consistent Rpp
        T = -2*Rpp/R - 2*Rp**2/R**2 - 4/R**2
        div_T = -2*Rpp/R - 2*Rp**2/R**2 - 1/R**2
        
        # Scalar equation: □φ - V'(φ) - ξφT + 2ξ∇_μ(φT^μ) = 0
        # □φ = (1/R²)∂_u(R²N∂_uφ) = (1/R²)(2R'Nφ' + R²N'φ' + R²Nφ'')
        dalembert_phi = (1/R**2) * (2*Rp*N*phip + R**2*0*phip + R**2*N*0)  # N' not yet computed
        dalembert_phi = (1/R**2) * (2*Rp*N*phip + R**2*N*0)  # Simplified
        
        F_prime = 2 * self.xi * phi
        
        # φ'' from scalar equation (simplified, need full expression)
        phipp = V_prime + 0.5 * F_prime * T - self.xi * phi * div_T  # Placeholder
        
        # N' from uu component or asymptotic condition
        # For wormhole: N → 1 at infinity, smooth transition at throat
        Np = 0.0  # Simplified - needs full uu equation
        
        return [Rp, Rpp, phip, phipp, Np]
    
    def solve_bvp(self, u_min=-10, u_max=10, n_points=500):
        """
        Solve RGS-EMS background using simplified ODE integration.
        
        Due to BVP complexity with algebraic constraints, use direct ODE integration
        with shooting to find boundary conditions.
        """
        try:
            # Simple ODE integration from initial conditions
            u = np.linspace(u_min, u_max, n_points)
            
            # Use Morris-Thorne-like profile as baseline (R(0)=2.31 from FORMULAZIONE)
            # R(u) = sqrt(u^2 + b0^2) with b0 = 2.31
            b0 = 2.31
            R = np.sqrt(u**2 + b0**2)
            Rp = u / R
            Rpp = b0**2 / (u**2 + b0**2)**1.5
            
            # Scalar field: Gaussian centered at throat
            phi = 0.8 * np.exp(-u**2)
            phip = -2 * 0.8 * u * np.exp(-u**2)
            
            # Lapse function: N≈1 (flat redshift approximation)
            N = np.ones_like(u)
            Np = np.zeros_like(u)
            
            # Compute torsion quantities
            T = -2*Rpp/R - 2*Rp**2/R**2 - 4/R**2
            div_T = -2*Rpp/R - 2*Rp**2/R**2 - 1/R**2
            
            # Verify teleparallel identity
            identity_residual = self.geometry.verify_teleparallel_identity(R, u)
            
            # Scalar equation residual
            scalar_residual = self.scalar_field_equation(R, phi, u)
            scalar_L2 = np.sqrt(np.trapz(np.abs(scalar_residual)**2, u))
            
            # Hamiltonian constraint
            G_tt = self.geometry.G_tt(N, R, u)
            hamiltonian_residual = np.sqrt(np.trapz(np.abs(G_tt)**2, u))
            
            # Physical quantities
            m_MS = self.quantities.misner_sharp_mass(R, self.G)
            E_em, rho_em = self.quantities.electromagnetic_energy_regularized(R, self.Q, self.ell0, self.G)
            
            # ADM extraction
            adm_fit = self.adm.asymptotic_fit(R, m_MS)
            
            # Extract throat values (at u=0, which is the center of the array)
            center_idx = len(u) // 2
            
            return {
                "status": "success",
                "method": "analytic_profile_with_residuals",
                "u": u,
                "R": R,
                "phi": phi,
                "N": N,
                "m_MS": m_MS,
                "E_EM": E_em,
                "rho_EM": rho_em,
                "residuals": {
                    "teleparallel_identity": identity_residual,
                    "scalar_L2": scalar_L2,
                    "hamiltonian_L2": hamiltonian_residual
                },
                "ADM_fit": adm_fit,
                "boundary_values": {
                    "R_0": R[center_idx],  # Throat radius at u=0
                    "phi_0": phi[center_idx],  # Scalar field at throat
                    "N_0": N[center_idx],  # Lapse at throat
                    "R_prime_0": Rp[center_idx],  # Should be ~0 at throat
                    "phi_prime_0": phip[center_idx]  # Should be ~0 at throat
                },
                "note": "Using analytic Morris-Thorne-like profile with R(0)=2.31 from FORMULAZIONE BVP result. Full BVP requires teleparallel gravitational equation completion."
            }
            
        except Exception as e:
            return {"status": "error", "message": str(e)}


# ==============================================================================
# SECTION 5: HAMILTONIAN CONSTRAINT VERIFICATION
# ==============================================================================

class HamiltonianConstraint:
    """
    Independent verification of Hamiltonian constraint.
    
    The constraint must be verified INDEPENDENTLY from evolution equations,
    not implicitly imposed through profile interpolation.
    """
    
    @staticmethod
    def compute(N, R, u, sources=None):
        """
        Compute Hamiltonian constraint H = 0.
        
        For spherically symmetric static metric:
        H = G^t_t - 8πG * T^t_t = 0
        
        Args:
            N, R: metric functions
            u: coordinate grid
            sources: matter source T^t_t (if available)
        
        Returns:
            H(u) constraint violation
        """
        geometry = GeometricIdentities()
        G_tt = geometry.G_tt(N, R, u)
        
        if sources is None:
            # For vacuum case, constraint is G_tt = 0
            return G_tt
        else:
            # With matter: G_tt - 8πG * T^t_t = 0
            return G_tt - 8 * np.pi * sources
    
    @staticmethod
    def verify_independent(solution, source_term=None):
        """
        Verify that Hamiltonian constraint is satisfied independently.
        
        Args:
            solution: dict with N, R, u from solver
            source_term: T^t_t if available
        
        Returns:
            verification report
        """
        H = HamiltonianConstraint.compute(
            solution["N"], solution["R"], solution["u"], source_term
        )
        
        residual = np.sqrt(np.trapz(np.abs(H)**2, solution["u"]))
        max_violation = np.max(np.abs(H))
        
        report = {
            "L2_residual": residual,
            "max_violation": max_violation,
            "satisfied": residual < 1e-6,  # tolerance
            "H_profile": H
        }
        
        return report


# ==============================================================================
# SECTION 6: GAP DOCUMENTATION
# ==============================================================================

def document_gaps():
    """
    Document the gaps in the RGS-PSU action specification.

    CRITICAL FINDING:
    The source documents containing the actual RGS-PSU field equations are NOT present
    in the recovered archive. Specifically:
    - FORMULAZIONE.pdf (Percorso A 4D background with R(0)=0.8292)
    - siiiiii.pdf (Profile with R(0)=2)

    Without these documents, the field equations cannot be derived from first principles.
    """
    gaps = {
        "source_documents": {
            "status": "MISSING",
            "description": "Source documents with RGS-PSU field equations not found in archive",
            "missing_documents": [
                "FORMULAZIONE.pdf - Contains Percorso A 4D background (R(0)=0.8292, N(0)=1.4529, φ(0)=0.3450)",
                "siiiiii.pdf - Contains profile with R(0)=2, N(0)=0.9821, φ(0)=0.3450"
            ],
            "searched_locations": [
                "/home/fabianleo/Scaricati/",
                "/home/fabianleo/Scaricati/PRISM_KristAll_recovered_sources/"
            ],
            "action_required": "Locate these documents or reconstruct field equations from source"
        },
        "action_specification": {
            "status": "INCOMPLETE",
            "description": "RGS-PSU action not uniquely specified in recovered sources",
            "missing": [
                "F_RGS(T, K_i) function in teleparallel action",
                "Complete NLED Lagrangian and variational sources",
                "Coupling between teleparallel theory and EMS model",
                "Specific form of regularization in electromagnetic sector"
            ]
        },
        "background_separation": {
            "status": "THREE_SEPARATE_CONFIGURATIONS",
            "description": "Three distinct background configurations found - do not mix",
            "configurations": {
                "placeholder_solver": {
                    "R_0": 1.0,
                    "source": "Morris-Thorne-like placeholder in current script",
                    "status": "NOT_VALIDATED",
                    "note": "This is NOT an RGS-PSU solution"
                },
                "percorso_A_4D": {
                    "R_0": 0.8292,
                    "N_0": 1.4529,
                    "phi_0": 0.3450,
                    "source": "FORMULAZIONE.pdf (missing)",
                    "status": "CANNOT_SOLVE",
                    "note": "Source document not available"
                },
                "siiiiii_profile": {
                    "R_0": 2.0,
                    "N_0": 0.9821,
                    "phi_0": 0.3450,
                    "source": "siiiiii.pdf (missing)",
                    "status": "CANNOT_SOLVE",
                    "note": "Source document not available"
                }
            }
        },
        "field_equations": {
            "status": "NOT_DERIVABLE",
            "description": "Field equations cannot be derived without source documents",
            "current_implementation": "PLACEHOLDER - uses Morris-Thorne profile for numerical stability only",
            "required_action": "When source documents are recovered, replace field_equations() method with actual RGS-PSU equations"
        }
    }

    return gaps


# ==============================================================================
# SECTION 7: MAIN EXECUTION
# ==============================================================================

def main():
    """
    Main execution: RGS-EMS corrected analysis and article generation.
    """
    print("=" * 80)
    print("RGS-PSU CLOSURE - CORRECTED ANALYSIS")
    print("=" * 80)
    print()
    print("FOUR CRITICAL INCONSISTENCIES IDENTIFIED:")
    print("A. Teleparallel identity: R_LC + T - 2∇_μT^μ = -2R''/R (not zero)")
    print("B. Tetrad: e^a_μ = diag(N,1,R,Rsinθ) not diag(1,1,R,Rsinθ)")
    print("C. EM field: E(r)~Q/r⁴ does NOT recover Coulomb Q/r²")
    print("D. Misner-Sharp: m_MS(0)=R(0)/2=1.155, not 5.124")
    print()
    
    # Initialize solver
    solver = BackgroundSolver(G=1.0, Q=1.0, ell0=1.0, xi=0.5, V0=0.1, lam=1.0, phi_v=0.5)
    
    # Test with Morris-Thorne profile
    print("SECTION 1: CORRECTED TELEPARALLEL IDENTITY TEST")
    print("-" * 80)
    u_grid = np.linspace(-10, 10, 1000)
    b0 = 2.31
    R_test = np.sqrt(u_grid**2 + b0**2)
    
    identity_result = solver.geometry.verify_teleparallel_identity(R_test, u_grid)
    print(f"  Identity value: R_LC + T - 2∇_μT^μ = {identity_result['mean_value']:.6e}")
    print(f"  Theoretical: -2R''/R = {np.mean(identity_result['theoretical_value']):.6e}")
    print(f"  Match: {np.allclose(identity_result['mean_value'], np.mean(identity_result['theoretical_value']))}")
    print(f"  Status: IDENTITY IS NOT ZERO (corrected)")
    print()
    
    # Corrected Misner-Sharp
    print("SECTION 2: CORRECTED MISNER-SHARP MASS")
    print("-" * 80)
    m_MS = solver.quantities.misner_sharp_mass(R_test, solver.G)
    center_idx = len(u_grid) // 2
    print(f"  R(0) = {R_test[center_idx]:.6f}")
    print(f"  R'(0) = {np.gradient(R_test)[center_idx]:.6e}")
    print(f"  m_MS(0) = {m_MS[center_idx]:.6f} (corrected: R(0)/2)")
    print(f"  Previous value 5.124 was INCORRECT")
    print()
    
    # Write complete article
    print("SECTION 3: WRITING COMPLETE ARTICLE")
    print("-" * 80)
    article = {
        "title": "RGS-PSU Closure: Mathematical Analysis and Inconsistency Documentation",
        "abstract": {
            "background": "This paper presents a mathematical analysis of the RGS-PSU (Relational Geometric-Structure with Phase-Shift) teleparallel model coupled to an Einstein-Maxwell-Scalar (EMS) system. Four critical inconsistencies are identified in the current implementation, preventing convergence to a validated solution.",
            "corrections": "The teleparallel geometric identity R_LC + T - 2∇_μT^μ is shown to equal -2R''/R, not zero in general. The tetrad for metric ds²=-N²dt²+du²+R²dΩ² must be diag(N,1,R,Rsinθ), not diag(1,1,R,Rsinθ). The electromagnetic field E(r)=Qr²/(r⁴+ℓ₀⁴)^(3/2)~Q/r⁴ does not recover the Coulomb limit Q/r². The Misner-Sharp mass at throat R(0)=2.31 with R'(0)=0 must be m_MS(0)=1.155, not 5.124.",
            "conclusion": "The RGS-PSU framework requires reconciliation of teleparallel conventions, derivation of a consistent NLED Lagrangian, and completion of the gravitational equation for N(u) dynamics before a convergent BVP solution can be obtained."
        },
        "sections": [
            {
                "number": 1,
                "title": "Introduction",
                "content": [
                    "Overview of RGS-PSU program and three background configurations",
                    "R(0)=2.31 (Percorso A), R(0)=0.8292, R(0)=2 - must remain distinct",
                    "Teleparallel geometry vs Levi-Civita formulation",
                    "Motivation for regularization in electromagnetic sector"
                ]
            },
            {
                "number": 2,
                "title": "Four Critical Inconsistencies",
                "subsection": [
                    {
                        "title": "2.1 Teleparallel Identity Failure",
                        "content": "For ultrastatic metric ds²=-dt²+du²+R²dΩ², R_LC=-4R''/R+2(1-R'²)/R². With T=-2R''/R-2R'²/R²-4/R² and ∇_μT^μ=-2R''/R-2R'²/R²-1/R², the identity R_LC+T-2∇_μT^μ=-2R''/R≠0 in general. Only vanishes for special R(u)."
                    },
                    {
                        "title": "2.2 Tetrad Incompatibility",
                        "content": "For metric ds²=-N(u)²dt²+du²+R(u)²dΩ², the correct diagonal tetrad is e^a_μ=diag(N(u),1,R(u),R(u)sinθ). The implemented diag(1,1,R,Rsinθ) corresponds to N=1 and cannot describe general N(u) dynamics."
                    },
                    {
                        "title": "2.3 Electromagnetic Field Inconsistency",
                        "content": "The density ρ_EM=Q²/(2√(r⁴+ℓ₀⁴)) asymptotically recovers Q²/(2r²) as expected. However, the reported field E(r)=Qr²/(r⁴+ℓ₀⁴)^(3/2) behaves as Q/r⁴, not Q/r². These cannot both follow from the same NLED Lagrangian without explicit derivation."
                    },
                    {
                        "title": "2.4 Misner-Sharp Mass Discrepancy",
                        "content": "For R(0)=2.31 and R'(0)=0, the Misner-Sharp mass formula m_MS=R(1-R'²)/2 gives m_MS(0)=1.155. The reported value 5.124 came from an inconsistent implementation and is incorrect."
                    }
                ]
            },
            {
                "number": 3,
                "title": "Required Corrections",
                "content": [
                    "1. Establish teleparallel conventions: tetrad, spin connection, torsion definition",
                    "2. Derive field equations from corrected action, not from incompatible formulas",
                    "3. Derive NLED Lagrangian that yields consistent ρ_EM and E(r) with correct limits",
                    "4. Implement correct tetrad with N(u) dependence in all torsion calculations",
                    "5. Solve coupled BVP for R(u), φ(u), N(u) with corrected equations"
                ]
            },
            {
                "number": 4,
                "title": "Current Status",
                "content": [
                    "Action structure specified but conventions not reconciled",
                    "Scalar potential V(φ)=V₀+(λ/4)(φ²-φ_v²)² available",
                    "Teleparallel identity: INCOMPATIBLE (corrected to -2R''/R)",
                    "Electromagnetic sector: DENSITY REGULAR, FIELD INCORRECT",
                    "Profile R(0)=2.31: ANALYTIC, NOT BVP SOLUTION",
                    "Gravitational BVP: INCOMPLETE (N(u) equation missing)"
                ]
            },
            {
                "number": 5,
                "title": "Conclusion",
                "content": [
                    "RGS-PSU framework identified with four mathematical inconsistencies",
                    "Current implementation cannot be presented as convergent solution",
                    "Path forward requires teleparallel convention reconciliation and NLED derivation",
                    "Three R(0) configurations remain distinct pending mathematical equivalence proof"
                ]
            }
        ],
        "corrected_formulas": [
            "R_LC = -4R''/R + 2(1-R'²)/R² (ultrastatic metric)",
            "T = -2R''/R - 2R'²/R² - 4/R²",
            "∇_μT^μ = -2R''/R - 2R'²/R² - 1/R²",
            "R_LC + T - 2∇_μT^μ = -2R''/R (NOT zero)",
            "e^a_μ = diag(N(u), 1, R(u), R(u)sinθ) (correct tetrad)",
            "m_MS(0) = R(0)/2 when R'(0)=0",
            "ρ_EM(r) = Q²/(2√(r⁴+ℓ₀⁴)) → Q²/(2r²) (OK)",
            "E(r) must be derived from NLED, E(r)~Q/r² required"
        ],
        "open_problems": [
            "Derive teleparallel gravitational equation for N(u) dynamics",
            "Derive NLED Lagrangian yielding consistent ρ_EM and E(r)",
            "Reconcile teleparallel conventions with source documents",
            "Prove or disprove equivalence of three R(0) configurations",
            "Complete BVP solver for corrected coupled system"
        ]
    }
    
    # Save corrected results
    output = {
        "timestamp": "2026-10-28",
        "status": "CORRECTED",
        "inconsistencies_documented": 4,
        "identity_check": identity_result,
        "m_MS_corrected": float(m_MS[center_idx]),
        "article": article
    }
    
    with open("RGS_PSU_closure_results.json", "w") as f:
        json.dump(output, f, indent=2)
    
    with open("RGS_PSU_article_complete.json", "w") as f:
        json.dump(article, f, indent=2)
    
    print("  Results: RGS_PSU_closure_results.json")
    print("  Article: RGS_PSU_article_complete.json")
    print()
    
    print("=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print("✓ Four inconsistencies identified and documented")
    print("✓ Teleparallel identity corrected to -2R''/R")
    print("✓ Misner-Sharp mass corrected to 1.155")
    print("✓ Tetrad and EM field inconsistencies documented")
    print("✓ Complete article with corrections written")
    print("→ Implementation cannot claim convergence without these fixes")
    print("=" * 80)
    print("=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print("✓ Geometric inconsistency report generated")
    print("→ No converged coupled R(u), phi(u), N(u) BVP is claimed")
    print("→ Teleparallel conventions and field equations remain OPEN")
    print("=" * 80)


if __name__ == "__main__":
    main()
