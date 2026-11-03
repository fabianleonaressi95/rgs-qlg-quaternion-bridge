# ============================================================================
# 64B point charges + uniform neutralizing CUBIC background: robust Kaggle code
# ============================================================================
# PURPOSE
# -------
# Test, without tuning kappa or any target parameter, whether the fixed 4x4x4
# grid of 64 like charges is an equilibrium/stable configuration in a uniform
# background charge. The background density is fixed ONLY by charge neutrality.
#
# The literal geometry puts 56 of 64 sites on the background boundary. Therefore
# derivatives at eta=0 require a boundary convention. This script uses the
# well-posed regulated family L = L_sites + eta, eta>0, with
# rho_b(eta) = -N/(2 L)^3. The meaningful question is the eta -> 0+ limit.
# eta is NOT fitted: it is only a boundary regulator and must be extrapolated.
#
# UNITS
# -----
# a = 1, q = 1, and 1/(4*pi*epsilon0) = 1.
# Restore: E_phys = E_code * q^2/(4*pi*epsilon0*a)
#          H_phys = H_code * q^2/(4*pi*epsilon0*a^3)
#
# WHAT IT COMPUTES
# ----------------
# Eqq: point-point energy
# Eqb: point-background energy
# Ebb: background-background self-energy
# F: forces on all sites
# H: 3N x 3N Hessian in particle coordinates, H=Hqq+Hqb, Hbb=0
# Spectrum: N-, N0, N+ with an explicit tolerance
# Delta H: verifies whether it is proportional to I
# A1g breathing vector: Rayleigh quotient and eigenvector residual
# n=3,4,5 energy/scaling runs plus n=4 Hessian runs by default.
#
# REQUIREMENTS: numpy, scipy (standard Kaggle Python image)
# ============================================================================

import itertools
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import eigh

np.set_printoptions(precision=14, suppress=False, linewidth=180)

# -----------------------------------------------------------------------------
# User settings
# -----------------------------------------------------------------------------
QUAD_ORDER = 32         # Increase 32 -> 40 -> 48 for a convergence check.
HESS_STEP = 2.0e-3      # Central finite-difference step for Hess Phi.
BB_ORDER = 48           # Duffy quadrature order for background self-energy.
ETAS = [1.0e-1, 3.0e-2, 1.0e-2, 3.0e-3]
RUN_N4_HESSIAN = True   # n=4 requires 192x192 diagonalization.
RUN_N3_N5_HESSIAN = False  # Set True for the full robustness test (slower).
SAVE_JSON = True

# -----------------------------------------------------------------------------
# Geometry
# -----------------------------------------------------------------------------
def make_sites(n):
    """n^3 sites: (i-(n-1)/2,j-(n-1)/2,k-(n-1)/2), in units a=1."""
    u = np.arange(n, dtype=float) - 0.5 * (n - 1)
    return np.array(np.meshgrid(u, u, u, indexing="ij")).reshape(3, -1).T

# -----------------------------------------------------------------------------
# Gauss-Legendre nodes on [0,1]
# -----------------------------------------------------------------------------
def gl01(order):
    z, w = leggauss(order)
    return 0.5 * (z + 1.0), 0.5 * w

# -----------------------------------------------------------------------------
# Duffy integration for a rectangular box having the observation point as corner
# -----------------------------------------------------------------------------
def positive_octant_box_integrals(dx, dy, dz, order=32):
    r"""
    For dx,dy,dz > 0 calculate:

      I0 = integral_[0,dx]x[0,dy]x[0,dz] d^3u / |u|
      G  = integral_[0,dx]x[0,dy]x[0,dz] u d^3u / |u|^3

    The six Duffy sectors remove the 1/r singularity analytically.

    IMPORTANT: if one dimension is zero, the corresponding 3D volume is zero.
    At eta>0, all eight boxes for all charge sites have positive dimensions.
    """
    d = np.array([dx, dy, dz], dtype=float)
    if np.any(d <= 0.0):
        return 0.0, np.zeros(3)

    s_nodes, s_weights = gl01(order)
    S, T, U = np.meshgrid(s_nodes, s_nodes, s_nodes, indexing="ij")
    W = (
        s_weights[:, None, None]
        * s_weights[None, :, None]
        * s_weights[None, None, :]
    )

    vol = float(np.prod(d))
    I0 = 0.0
    G = np.zeros(3)

    # Partition according to the ordering of the normalized rectangular coords.
    # z_p0=s, z_p1=s*t, z_p2=s*u. Physical vector = s*A(t,u).
    for perm in itertools.permutations((0, 1, 2)):
        A = np.zeros((3,) + S.shape, dtype=float)
        A[perm[0]] = d[perm[0]]
        A[perm[1]] = d[perm[1]] * T
        A[perm[2]] = d[perm[2]] * U

        normA2 = np.sum(A * A, axis=0)
        normA = np.sqrt(normA2)

        # Gauss nodes are strictly inside (0,1), so normA should never be zero.
        # The masks guarantee no RuntimeWarning if roundoff/degenerate input occurs.
        invA = np.divide(
            1.0, normA, out=np.zeros_like(normA), where=normA > 1.0e-300
        )
        invA3 = invA * invA * invA

        # Jacobian = vol*s^2; |u|=s*|A| => contribution vol*s/|A|.
        I0 += np.sum(W * vol * S * invA)

        # Jacobian*u_i/|u|^3 = vol*A_i/|A|^3.
        for i in range(3):
            G[i] += np.sum(W * vol * A[i] * invA3)

    return float(I0), G

# -----------------------------------------------------------------------------
# Cube potential and gradient
# -----------------------------------------------------------------------------
def cube_phi_grad(x, L, order=32):
    r"""
    Compute for a uniformly filled cube [-L,L]^3:

      Phi(x) = integral_cube d^3r / |r-x|
      grad_x Phi(x) = integral_cube (r-x)/|r-x|^3 d^3r

    x must lie strictly inside the cube. The caller enforces this through eta>0.
    """
    x = np.asarray(x, dtype=float)
    if np.any(np.abs(x) >= L):
        raise ValueError(
            "Observation point touches/exits the cube. Use eta>0 so all sites are interior."
        )

    phi = 0.0
    grad = np.zeros(3)

    # Split the cube into eight rectangular boxes which have x as a common corner.
    for sx in (-1.0, 1.0):
        dx = (x[0] + L) if sx < 0.0 else (L - x[0])
        for sy in (-1.0, 1.0):
            dy = (x[1] + L) if sy < 0.0 else (L - x[1])
            for sz in (-1.0, 1.0):
                dz = (x[2] + L) if sz < 0.0 else (L - x[2])

                val, g_positive = positive_octant_box_integrals(dx, dy, dz, order)
                phi += val

                # Source displacement r-x has signs (sx,sy,sz).
                grad += np.array([sx, sy, sz]) * g_positive

    return float(phi), grad

# -----------------------------------------------------------------------------
# Background potential Hessian: finite difference of the regularized gradient
# -----------------------------------------------------------------------------
def cube_hessian_phi(x, L, order=32, h=2.0e-3):
    """Hessian of Phi at a strictly interior point using central differences."""
    x = np.asarray(x, dtype=float)
    if np.any(np.abs(x) + h >= L):
        raise ValueError(
            "Finite-difference stencil crosses the background boundary. Increase eta or reduce h."
        )

    H = np.zeros((3, 3))
    basis = np.eye(3)

    for j in range(3):
        _, gp = cube_phi_grad(x + h * basis[j], L, order)
        _, gm = cube_phi_grad(x - h * basis[j], L, order)
        H[:, j] = (gp - gm) / (2.0 * h)

    return 0.5 * (H + H.T)

# -----------------------------------------------------------------------------
# Point-point Coulomb energy, force, Hessian
# -----------------------------------------------------------------------------
def qq_energy_force_hessian(X):
    """Units q^2/(4*pi*eps0)=1."""
    N = len(X)
    E = 0.0
    F = np.zeros((N, 3))
    H = np.zeros((3 * N, 3 * N))
    I3 = np.eye(3)

    for A in range(N):
        for B in range(A + 1, N):
            d = X[A] - X[B]
            r2 = float(d @ d)
            r = np.sqrt(r2)

            E += 1.0 / r

            # Physical force on A from B.
            f = d / r**3
            F[A] += f
            F[B] -= f

            # Hessian contribution for pair energy 1/r.
            K = 3.0 * np.outer(d, d) / r**5 - I3 / r**3
            a = slice(3 * A, 3 * A + 3)
            b = slice(3 * B, 3 * B + 3)
            H[a, a] += K
            H[b, b] += K
            H[a, b] -= K
            H[b, a] -= K

    return float(E), F, H

# -----------------------------------------------------------------------------
# Background-background energy: Duffy quadrature after convolution reduction
# -----------------------------------------------------------------------------
def background_convolution_integral(L, order=48):
    r"""
    J(L) = int_V dr int_V dr' 1/|r-r'|,
    V=[-L,L]^3.

    By convolution:
      J = integral_[−2L,2L]^3 prod_i(2L-|u_i|)/|u| d^3u.

    Symmetry gives eight positive octants. Duffy sectors regularize u=0.
    """
    nodes, weights = gl01(order)
    S, T, U = np.meshgrid(nodes, nodes, nodes, indexing="ij")
    W = weights[:, None, None] * weights[None, :, None] * weights[None, None, :]

    total_positive_octant = 0.0

    for perm in itertools.permutations((0, 1, 2)):
        z = np.zeros((3,) + S.shape)
        z[perm[0]] = S
        z[perm[1]] = S * T
        z[perm[2]] = S * U

        A = np.zeros((3,) + S.shape)
        A[perm[0]] = 1.0
        A[perm[1]] = T
        A[perm[2]] = U
        normA = np.sqrt(np.sum(A * A, axis=0))
        invA = np.divide(1.0, normA, out=np.zeros_like(normA), where=normA > 1.0e-300)

        # u_i=2L*z_i; Jacobian=(2L)^3*S^2;
        # 1/|u|=1/(2L*S*|A|);
        # convolution weight=(2L)^3 prod(1-z_i).
        integrand = (2.0 * L)**5 * S * np.prod(1.0 - z, axis=0) * invA
        total_positive_octant += np.sum(W * integrand)

    return float(8.0 * total_positive_octant)

# -----------------------------------------------------------------------------
# Build the full model for a given n and eta
# -----------------------------------------------------------------------------
def solve_model(n, eta, quad_order=32, hess_step=2e-3, bb_order=48, do_hessian=True):
    X = make_sites(n)
    N = len(X)
    L_sites = 0.5 * (n - 1)
    L = L_sites + float(eta)

    if eta <= hess_step:
        raise ValueError("Need eta > HESS_STEP so central differences remain inside the cube.")

    rho_b = -N / (2.0 * L)**3

    Eqq, Fqq, Hqq = qq_energy_force_hessian(X)

    phis = np.empty(N)
    grad_phi = np.empty((N, 3))

    for A, x in enumerate(X):
        phis[A], grad_phi[A] = cube_phi_grad(x, L, quad_order)

    Eqb = rho_b * float(np.sum(phis))
    Fqb = -rho_b * grad_phi

    Jbb = background_convolution_integral(L, bb_order)
    Ebb = 0.5 * rho_b**2 * Jbb

    result = {
        "n": int(n),
        "N": int(N),
        "eta": float(eta),
        "L_sites": float(L_sites),
        "L_background": float(L),
        "rho_b": float(rho_b),
        "Eqq": float(Eqq),
        "Eqb": float(Eqb),
        "Ebb": float(Ebb),
        "Etot": float(Eqq + Eqb + Ebb),
        "Fqq_max": float(np.max(np.linalg.norm(Fqq, axis=1))),
        "Fqb_max": float(np.max(np.linalg.norm(Fqb, axis=1))),
        "Ftotal_max": float(np.max(np.linalg.norm(Fqq + Fqb, axis=1))),
        "Ftotal_rms_component": float(np.sqrt(np.mean((Fqq + Fqb)**2))),
        "Jbb": float(Jbb),
        "w_classical": 1.0 / 3.0,
    }

    if not do_hessian:
        return result

    Hqb = np.zeros_like(Hqq)
    for A, x in enumerate(X):
        block = cube_hessian_phi(x, L, quad_order, hess_step)
        Hqb[3 * A:3 * A + 3, 3 * A:3 * A + 3] = rho_b * block

    # Ebb has no dependence on particle positions.
    Htot = Hqq + Hqb
    eigvals = eigh(Htot, eigvals_only=True, check_finite=False)
    tol = max(1.0e-8, 1.0e-7 * float(np.max(np.abs(eigvals))))

    # The background is fixed in space: translations of particles alone are NOT
    # exact zero modes. These Rayleigh quotients document this rather than assume it.
    translation_rayleigh = []
    for axis in range(3):
        v = np.zeros(3 * N)
        v[axis::3] = 1.0
        translation_rayleigh.append(float(v @ Htot @ v / (v @ v)))

    # Candidate global breathing vector: delta x_A proportional to x_A.
    breathing = X.reshape(-1).copy()
    breathing /= np.linalg.norm(breathing)
    lambda_breath = float(breathing @ Htot @ breathing)
    breathing_residual = float(np.linalg.norm(Htot @ breathing - lambda_breath * breathing))

    mean_shift = float(np.trace(Hqb) / (3 * N))
    identity_residual = float(np.linalg.norm(Hqb - mean_shift * np.eye(3 * N)))

    result.update({
        "eig_min": float(eigvals[0]),
        "eig_max": float(eigvals[-1]),
        "eig_tol": float(tol),
        "Nneg": int(np.sum(eigvals < -tol)),
        "Nzero": int(np.sum(np.abs(eigvals) <= tol)),
        "Npos": int(np.sum(eigvals > tol)),
        "translation_rayleigh": translation_rayleigh,
        "breathing_rayleigh": lambda_breath,
        "breathing_residual": breathing_residual,
        "deltaH_frobenius": float(np.linalg.norm(Hqb)),
        "deltaH_identity_shift": mean_shift,
        "deltaH_identity_residual": identity_residual,
        "eigenvalues": eigvals.tolist(),
    })

    return result

# -----------------------------------------------------------------------------
# Basic tests of potential/gradient symmetry
# -----------------------------------------------------------------------------
def potential_sanity_checks():
    print("\n=== Potential/gradient sanity checks ===")
    L = 1.6
    points = [
        np.array([0.0, 0.0, 0.0]),
        np.array([0.5, 0.5, 0.5]),
        np.array([1.2, 0.4, 0.2]),
    ]
    for x in points:
        phi, g = cube_phi_grad(x, L, QUAD_ORDER)
        _, gm = cube_phi_grad(-x, L, QUAD_ORDER)
        print("x =", x)
        print("Phi =", phi)
        print("grad Phi =", g)
        print("grad Phi(x)+grad Phi(-x) norm =", np.linalg.norm(g + gm))
        print("finite?", np.isfinite(phi) and np.all(np.isfinite(g)))

# -----------------------------------------------------------------------------
# Convergence checks
# -----------------------------------------------------------------------------
def convergence_checks():
    print("\n=== Background self-energy convergence: n=4, eta=0.03 ===")
    n = 4
    N = n**3
    L = 1.5 + 0.03
    rho = -N / (2.0 * L)**3
    previous = None
    for order in (16, 24, 32, 40, 48, 64):
        J = background_convolution_integral(L, order)
        Ebb = 0.5 * rho**2 * J
        delta = None if previous is None else abs(Ebb - previous)
        print(f"order={order:3d}   Ebb={Ebb:.14f}   delta={delta}")
        previous = Ebb

    print("\n=== Gradient quadrature convergence: n=4 corner site, eta=0.03 ===")
    x = np.array([1.5, 1.5, 1.5])
    for order in (16, 24, 32, 40, 48):
        phi, g = cube_phi_grad(x, L, order)
        print(f"order={order:3d}   Phi={phi:.14f}   grad={g}")

# -----------------------------------------------------------------------------
# Driver
# -----------------------------------------------------------------------------
def compact_print(r):
    keys = [
        "n", "N", "eta", "L_background", "rho_b",
        "Eqq", "Eqb", "Ebb", "Etot",
        "Fqq_max", "Fqb_max", "Ftotal_max", "Ftotal_rms_component",
        "w_classical",
    ]
    for k in keys:
        print(f"{k}: {r[k]}")
    if "eig_min" in r:
        print(f"eig_min: {r['eig_min']}")
        print(f"eig_max: {r['eig_max']}")
        print(f"Nneg, Nzero, Npos: {r['Nneg']}, {r['Nzero']}, {r['Npos']}")
        print(f"eig_tol: {r['eig_tol']}")
        print(f"translation_rayleigh: {r['translation_rayleigh']}")
        print(f"breathing_rayleigh: {r['breathing_rayleigh']}")
        print(f"breathing_residual: {r['breathing_residual']}")
        print(f"deltaH_frobenius: {r['deltaH_frobenius']}")
        print(f"deltaH_identity_shift: {r['deltaH_identity_shift']}")
        print(f"deltaH_identity_residual: {r['deltaH_identity_residual']}")


if __name__ == "__main__":
    potential_sanity_checks()
    convergence_checks()

    all_results = []

    # Energy/force robustness for 3^3, 4^3, 5^3.
    for n in (3, 4, 5):
        print("\n" + "=" * 78)
        print(f"ENERGY + FORCE TEST: n={n}, N={n**3}")
        print("=" * 78)
        for eta in ETAS:
            # Skip too-small eta for central Hessian stencil only; force remains well-defined.
            do_h = (n == 4 and RUN_N4_HESSIAN and eta > HESS_STEP)
            if n == 3 and RUN_N3_N5_HESSIAN and eta > HESS_STEP:
                do_h = True
            if n == 5 and RUN_N3_N5_HESSIAN and eta > HESS_STEP:
                do_h = True

            print(f"\n--- n={n}, eta={eta} ---")
            r = solve_model(
                n=n,
                eta=eta,
                quad_order=QUAD_ORDER,
                hess_step=HESS_STEP,
                bb_order=BB_ORDER,
                do_hessian=do_h,
            )
            compact_print(r)
            all_results.append(r)

    # Exact scaling statement, independent of numerical constants.
    print("\n" + "=" * 78)
    print("ANALYTIC SCALING")
    print("=" * 78)
    print("At fixed n and fixed dimensionless eta:")
    print("Eqq(a), Eqb(a), Ebb(a) each scale as a^(-1).")
    print("Thus Etot(a)=K/a, V(a)=[2(L_sites+eta)a]^3, and p=rho_energy/3.")
    print("Hence w_classical=+1/3 exactly; this does not imply mechanical equilibrium.")
    print("If a stable spectrum existed and m were an a-independent external input:")
    print("H(a)~a^(-3), omega~a^(-3/2), EZPE~a^(-3/2), so w_ZPE=+1/2.")
    print("hbar and m remain INPUTS; do not compute a ZPE if Nneg>0.")

    if SAVE_JSON:
        clean = []
        for r in all_results:
            # Keep JSON small: eigenvalues are available in live output but omitted from summary.
            q = {k: v for k, v in r.items() if k != "eigenvalues"}
            clean.append(q)
        out = Path("64B_cubic_jellium_results.json")
        out.write_text(json.dumps(clean, indent=2))
        print(f"\nSaved summary: {out}")
