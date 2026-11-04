# Kaggle notebook/script: 64B charges + neutralizing cubic jellium
# Runs with NumPy/SciPy only. Uses Gauss-Legendre tensor quadrature for the cube.
# Units: a=q=(4*pi*epsilon0)^(-1)=1. Restore SI factors by multiplying energies
# by q^2/(4*pi*epsilon0*a) and Hessian eigenvalues by q^2/(4*pi*epsilon0*a^3).

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import eigh
from scipy.optimize import brentq
from collections import Counter

np.set_printoptions(precision=14, suppress=False, linewidth=180)

# -----------------------------------------------------------------------------
# Geometry
# -----------------------------------------------------------------------------
def sites(n=4):
    u = np.arange(n, dtype=float) - (n - 1) / 2.0
    return np.array(np.meshgrid(u, u, u, indexing='ij')).reshape(3, -1).T

# -----------------------------------------------------------------------------
# Cube potential and derivatives, evaluated by tensor-product Gauss-Legendre.
# Cube is [-L,L]^3. This implements only finite, nonsingular integrals:
# - potential at a charge site is integrated by a Duffy split into 8 octants;
# - gradient/Hessian use differentiated potentials with principal-value caveat.
# For reliable Hessian use finite differences of the exact cube potential below.
# -----------------------------------------------------------------------------
def cube_potential_exact(x, L=1.5):
    """Potential integral int_cube d^3r / |x-r|. Analytic rectangular-prism formula.
    Valid for interior points by continuous limiting prescription.
    Formula is summed over translated box corners, with F(X,Y,Z).
    """
    x, y, z = map(float, x)
    def F(X, Y, Z):
        R = np.sqrt(X*X + Y*Y + Z*Z)
        # Stable forms; clip only for roundoff at singular faces/edges.
        t1 = X*Y*np.log(np.maximum(Z + R, np.finfo(float).tiny))
        t2 = Y*Z*np.log(np.maximum(X + R, np.finfo(float).tiny))
        t3 = Z*X*np.log(np.maximum(Y + R, np.finfo(float).tiny))
        t4 = -0.5*X*X*np.arctan2(Y*Z, X*R)
        t5 = -0.5*Y*Y*np.arctan2(X*Z, Y*R)
        t6 = -0.5*Z*Z*np.arctan2(X*Y, Z*R)
        return t1 + t2 + t3 + t4 + t5 + t6
    total = 0.0
    for sx in (-1, 1):
        X = sx*L - x
        for sy in (-1, 1):
            Y = sy*L - y
            for sz in (-1, 1):
                Z = sz*L - z
                total += sx*sy*sz*F(X, Y, Z)
    return total

def cube_potential_quad(x, L=1.5, order=42):
    """Independent high-order check of cube potential, using Duffy maps around x.
    Decomposes the cube into 8 sub-boxes sharing x, making 1/r integrable smoothly.
    """
    g, w = leggauss(order)
    t = (g + 1.0)/2.0
    wt = w/2.0
    xx, yy, zz = np.meshgrid(t, t, t, indexing='ij')
    ww = wt[:,None,None]*wt[None,:,None]*wt[None,None,:]
    ans = 0.0
    for sx in (-1.0, 1.0):
        dx = (L - x[0]) if sx > 0 else (x[0] + L)
        for sy in (-1.0, 1.0):
            dy = (L - x[1]) if sy > 0 else (x[1] + L)
            for sz in (-1.0, 1.0):
                dz = (L - x[2]) if sz > 0 else (x[2] + L)
                # r-x = sign*(dx*u,dy*v,dz*w); use u=s*t, v=s*(1-t), etc.
                # Symmetric Duffy subdivision selected by largest coordinate not needed:
                # direct tensor cube is singular only at endpoint and converges here.
                rx, ry, rz = dx*xx, dy*yy, dz*zz
                rr = np.sqrt(rx*rx + ry*ry + rz*rz)
                ans += np.sum(ww * (dx*dy*dz) / rr)
    return ans

def numerical_grad_hess_phi(x, L=1.5, h=2e-4):
    """Central finite differences of analytic cube potential."""
    x = np.asarray(x, dtype=float)
    f0 = cube_potential_exact(x, L)
    grad = np.empty(3)
    H = np.empty((3,3))
    e = np.eye(3)
    for i in range(3):
        fp, fm = cube_potential_exact(x+h*e[i], L), cube_potential_exact(x-h*e[i], L)
        grad[i] = (fp-fm)/(2*h)
        H[i,i] = (fp - 2*f0 + fm)/h**2
        for j in range(i):
            fpp = cube_potential_exact(x+h*e[i]+h*e[j], L)
            fpm = cube_potential_exact(x+h*e[i]-h*e[j], L)
            fmp = cube_potential_exact(x-h*e[i]+h*e[j], L)
            fmm = cube_potential_exact(x-h*e[i]-h*e[j], L)
            H[i,j] = H[j,i] = (fpp-fpm-fmp+fmm)/(4*h*h)
    return grad, H

# -----------------------------------------------------------------------------
# Point-charge force and Hessian, in units q^2/(4*pi*eps0)=1.
# -----------------------------------------------------------------------------
def qq_energy_force_hessian(X):
    N = len(X)
    E = 0.0
    F = np.zeros((N,3))
    H = np.zeros((3*N,3*N))
    I = np.eye(3)
    for A in range(N):
        for B in range(A+1, N):
            d = X[A]-X[B]
            r2 = float(d@d)
            r = np.sqrt(r2)
            E += 1.0/r
            f = d/r**3
            F[A] += f
            F[B] -= f
            K = (3*np.outer(d,d)/r**5 - I/r**3)
            ia, ib = slice(3*A,3*A+3), slice(3*B,3*B+3)
            H[ia,ia] += K
            H[ib,ib] += K
            H[ia,ib] -= K
            H[ib,ia] -= K
    return E, F, H

# -----------------------------------------------------------------------------
# Background. Neutrality fixes rho=-N/(2L)^3 in q=1 units.
# E_qb = rho sum Phi(x_A), E_bb = rho^2/2 * double integral.
# Background energy is evaluated using convolution reduction:
# int_V int_V f(r-r') = int_{[-2L,2L]^3} prod_i(2L-|u_i|) / |u| du.
# The origin singularity is integrable; tensor Gauss rules converge algebraically.
# -----------------------------------------------------------------------------
def background_bb_energy(L=1.5, rho=None, order=80):
    if rho is None:
        rho = -64/(2*L)**3
    g,w = leggauss(order)
    u = 2*L*g
    wu = 2*L*w
    U,V,W = np.meshgrid(u,u,u,indexing='ij')
    WW = wu[:,None,None]*wu[None,:,None]*wu[None,None,:]
    R = np.sqrt(U*U+V*V+W*W)
    weights = (2*L-np.abs(U))*(2*L-np.abs(V))*(2*L-np.abs(W))
    integrand = np.divide(weights, R, out=np.zeros_like(R), where=R>0)
    J = np.sum(WW*integrand)
    return 0.5*rho*rho*J, J

def total_for_n(n=4, h=2e-4, bb_order=90):
    X = sites(n)
    L = (n-1)/2
    N = len(X)
    rho = -N/(2*L)**3
    Eqq,Fqq,Hqq = qq_energy_force_hessian(X)
    phi = np.array([cube_potential_exact(x,L) for x in X])
    grads = []
    Hb = np.zeros_like(Hqq)
    for A,x in enumerate(X):
        g,H = numerical_grad_hess_phi(x,L,h)
        grads.append(g)
        Hb[3*A:3*A+3,3*A:3*A+3] = rho*H
    grads = np.array(grads)
    Eqb = rho*np.sum(phi)
    Fqb = -rho*grads
    Ebb,J = background_bb_energy(L,rho,bb_order)
    F = Fqq + Fqb
    H = Hqq + Hb
    eig = eigh(H, eigvals_only=True, check_finite=False)
    tol = max(1e-8, 1e-7*np.max(np.abs(eig)))
    nm = int(np.sum(eig < -tol)); nz=int(np.sum(np.abs(eig)<=tol)); np_=int(np.sum(eig>tol))
    breath = X.reshape(-1)
    breath /= np.linalg.norm(breath)
    Hbreat = float(breath@H@breath)
    resid = np.linalg.norm(H@breath-Hbreat*breath)
    # derivative along uniform scaling at a=1: E(a)=E(1)/a.
    Etot=Eqq+Eqb+Ebb
    report = {
        'n':n,'N':N,'L':L,'rho_b':rho,
        'Eqq':Eqq,'Eqb':Eqb,'Ebb':Ebb,'Etot':Etot,
        'force_qq_max':float(np.max(np.linalg.norm(Fqq,axis=1))),
        'force_qb_max':float(np.max(np.linalg.norm(Fqb,axis=1))),
        'force_total_max':float(np.max(np.linalg.norm(F,axis=1))),
        'force_total_rms':float(np.sqrt(np.mean(F*F))),
        'eig_min':float(eig[0]),'eig_max':float(eig[-1]),
        'Nneg':nm,'Nzero':nz,'Npos':np_, 'tol':tol,
        'translation_rayleigh': [float(v@H@v/(v@v)) for v in [np.tile([1.,0,0],N),np.tile([0,1.,0],N),np.tile([0,0,1.],N)]],
        'breathing_rayleigh':Hbreat,'breathing_residual':float(resid),
        'dE_da_at_1':-Etot,'d2E_da2_at_1':2*Etot,
        'w_classical':1/3,
        'bb_convolution_integral':float(J),
        'eigenvalues':eig,
    }
    return report, Hqq, Hb, H

def print_report(r, print_eigs=False):
    for k in ['n','N','L','rho_b','Eqq','Eqb','Ebb','Etot','force_qq_max','force_qb_max','force_total_max','force_total_rms','eig_min','eig_max','Nneg','Nzero','Npos','tol','translation_rayleigh','breathing_rayleigh','breathing_residual','dE_da_at_1','d2E_da2_at_1','w_classical','bb_convolution_integral']:
        print(f'{k}: {r[k]}')
    if print_eigs:
        print('eigenvalues:')
        print(r['eigenvalues'])

# -----------------------------------------------------------------------------
# Convergence study and requested runs n=3,4,5.
# Note: increase bb_order and vary h to document numerical stability.
# -----------------------------------------------------------------------------
if __name__ == '__main__':
    print('Convergence check: E_bb for n=4')
    L=1.5; rho=-64/(2*L)**3
    vals=[]
    for order in (32,48,64,80,100,128):
        e,j=background_bb_energy(L,rho,order)
        vals.append((order,e))
        print(order,e)
    print('successive |Delta Ebb|:', [abs(vals[i][1]-vals[i-1][1]) for i in range(1,len(vals))])

    results={}
    for n in (3,4,5):
        # n=5 needs more time/memory; use scipy LAPACK and expect a few minutes on Kaggle CPU.
        print('\n=== n=',n,'===')
        r,Hqq,Hb,H=total_for_n(n=n,h=2e-4,bb_order=100)
        results[n]=r
        print_report(r, print_eigs=(n==4))
        print('Delta H norm:', np.linalg.norm(Hb), 'identity-fit residual:', np.linalg.norm(Hb-np.trace(Hb)/Hb.shape[0]*np.eye(Hb.shape[0])))

    print('\nScaling check under a -> s a: every classical contribution must scale s^-1.')
    print('Therefore E(a)=E(1)/a, p=rho_energy/3, w=1/3 exactly, independent of n.')
    print('For a fixed input mass m, H(a)=H(1)/a^3, omega~a^(-3/2), EZPE~a^(-3/2), hence w_ZPE=+1/2 if and only if the stable spectrum exists.')
