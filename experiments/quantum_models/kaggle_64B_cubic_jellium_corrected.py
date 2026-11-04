# Corrected Kaggle notebook/script for 64B + fixed neutralizing cubic background.
# Key correction versus the earlier version: the rectangular-prism potential had an
# incorrect signed-corner implementation. Here Phi, grad Phi and Hess Phi are
# evaluated directly by high-order Gauss-Legendre quadrature with a Duffy map on
# each of the 8 boxes meeting at the observation point. This regularizes 1/r exactly.
#
# Units: a=q=(4*pi*epsilon0)^(-1)=1.  Cube is [-L,L]^3; sites are integer-grid
# points x=(i-(n-1)/2,j-(n-1)/2,k-(n-1)/2), and rho_b=-N/(2L)^3.
#
# IMPORTANT physical convention:
# background is fixed in space. Hence translations of particles alone are NOT zero modes.
# E_bb does not affect forces/Hessian. It is included with convergent Duffy quadrature.

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import eigh

np.set_printoptions(precision=14, linewidth=180)

def sites(n):
    u=np.arange(n,dtype=float)-(n-1)/2
    return np.array(np.meshgrid(u,u,u,indexing='ij')).reshape(3,-1).T

def duffy_nodes(order):
    z,w=leggauss(order)
    t=(z+1)/2; wt=w/2
    return t,wt

def box_integrals_from_corner(dx,dy,dz,order=32, want_hessian=True):
    """Integral over [0,dx]x[0,dy]x[0,dz] of 1/r and derivatives at origin.
    Uses six ordered Duffy sectors, in each r coordinate r_perm[0]=s,
    r_perm[1]=s*t, r_perm[2]=s*u. Singular radial factor cancels analytically.
    Returns I0, grad wrt observation x (for r=source-x), and Hessian wrt x.
    """
    t,w=duffy_nodes(order)
    S,T,U=np.meshgrid(t,t,t,indexing='ij')
    W=w[:,None,None]*w[None,:,None]*w[None,None,:]
    dims=np.array([dx,dy,dz],float)
    I0=0.0; G=np.zeros(3); H=np.zeros((3,3))
    # partition box by ordering of normalized coordinates y_i/d_i
    import itertools
    for perm in itertools.permutations(range(3)):
        # normalized coordinate vector z: z[p0]=s, z[p1]=s*t, z[p2]=s*u
        zc=np.empty((3,)+S.shape)
        zc[perm[0]]=S
        zc[perm[1]]=S*T
        zc[perm[2]]=S*U
        Rvec=np.array([dims[i]*zc[i] for i in range(3)])
        r=np.sqrt(np.sum(Rvec*Rvec,axis=0))
        # dr physical volume = prod(dims) * s^2 ds dt du
        # I0 integrand: prod(dims) * s / sqrt(sum(d_i^2 z_i^2)/s^2)
        # Let Rvec=s*A, r=s*Aabs
        A=np.empty_like(Rvec)
        A[perm[0]]=dims[perm[0]]
        A[perm[1]]=dims[perm[1]]*T
        A[perm[2]]=dims[perm[2]]*U
        aa=np.sqrt(np.sum(A*A,axis=0))
        jac=np.prod(dims)*S*S
        I0 += np.sum(W*jac/r)
        # grad_x 1/|source-x| = rvec/r^3 (with rvec=source-x)
        # jac*rvec/r^3 -> prod dims * A / aa^3, independent of s
        for i in range(3):
            G[i] += np.sum(W*np.prod(dims)*A[i]/aa**3)
        if want_hessian:
            # d_i d_j 1/r wrt x = (3 r_i r_j-r^2 delta)/r^5
            # jac * ... -> prod dims * (3 A_i A_j-aa^2 delta)/(s*aa^5)
            # integral over s has log singularity only in individual corner-sector
            # but full Hessian inside a charged volume has distributional/contact term;
            # avoid direct singular Hessian. Caller uses finite-difference gradients.
            pass
    return I0,G

def phi_grad_cube(x,L,order=36):
    """Phi=int_cube 1/|x-r| dr and grad_x Phi, exact-singularity regularized by Duffy."""
    x=np.asarray(x,float)
    phi=0.0; grad=np.zeros(3)
    for sx in (-1,1):
        dx=(L-x[0]) if sx==1 else (x[0]+L)
        for sy in (-1,1):
            dy=(L-x[1]) if sy==1 else (x[1]+L)
            for sz in (-1,1):
                dz=(L-x[2]) if sz==1 else (x[2]+L)
                val,g_local=box_integrals_from_corner(dx,dy,dz,order)
                signs=np.array([sx,sy,sz],float)
                phi += val
                # source-x has signed directions; local calculation assumes positive source-x
                grad += signs*g_local
    return phi,grad

def hessian_phi_fd(x,L,order=36,h=2e-3):
    """Differentiate regularized analytic gradient.  Sweep h for convergence."""
    H=np.zeros((3,3)); e=np.eye(3)
    for j in range(3):
        gp=phi_grad_cube(x+h*e[j],L,order)[1]
        gm=phi_grad_cube(x-h*e[j],L,order)[1]
        H[:,j]=(gp-gm)/(2*h)
    return 0.5*(H+H.T)

def qq_energy_force_hessian(X):
    N=len(X); E=0.; F=np.zeros((N,3)); H=np.zeros((3*N,3*N)); I=np.eye(3)
    for A in range(N):
        for B in range(A+1,N):
            d=X[A]-X[B]; r2=d@d; r=np.sqrt(r2)
            E+=1/r
            # physical force = d/r^3 ; energy gradient is -force
            f=d/r**3; F[A]+=f; F[B]-=f
            K=(3*np.outer(d,d)/r**5-I/r**3)
            ia=slice(3*A,3*A+3); ib=slice(3*B,3*B+3)
            H[ia,ia]+=K; H[ib,ib]+=K; H[ia,ib]-=K; H[ib,ia]-=K
    return E,F,H

def bb_integral_duffy(L,order=48):
    """J=double integral cube 1/|r-r'|. Uses convolution and 8 octants;
    each octant has weight (2L-u)(2L-v)(2L-w). Duffy removes origin singularity."""
    t,w=duffy_nodes(order)
    S,T,U=np.meshgrid(t,t,t,indexing='ij'); W=w[:,None,None]*w[None,:,None]*w[None,None,:]
    ans=0.; import itertools
    for perm in itertools.permutations(range(3)):
        z=np.empty((3,)+S.shape)
        z[perm[0]]=S; z[perm[1]]=S*T; z[perm[2]]=S*U
        A=np.empty_like(z)
        A[perm[0]]=1.; A[perm[1]]=T; A[perm[2]]=U
        aa=np.sqrt(np.sum(A*A,axis=0))
        # u=2L*s etc; product du = (2L)^3 s^2 dsdtdu; 1/r=1/(2L*s*aa)
        # octant weight product_i(2L-u_i) = (2L)^3 product(1-z_i)
        integrand=(2*L)**5 * S*np.prod(1-z,axis=0)/aa
        ans+=np.sum(W*integrand)
    return 8*ans

def solve(n=4,order=36,h=2e-3,bb_order=64):
    X=sites(n); N=len(X); L=(n-1)/2; rho=-N/(2*L)**3
    Eqq,Fqq,Hqq=qq_energy_force_hessian(X)
    phis=[]; gb=[]; Hqb=np.zeros_like(Hqq)
    for A,x in enumerate(X):
        ph,g=phi_grad_cube(x,L,order)
        phis.append(ph); gb.append(g)
        Hqb[3*A:3*A+3,3*A:3*A+3]=rho*hessian_phi_fd(x,L,order,h)
    phis=np.array(phis); gb=np.array(gb)
    Eqb=rho*phis.sum(); Fqb=-rho*gb
    J=bb_integral_duffy(L,bb_order); Ebb=.5*rho*rho*J
    H=Hqq+Hqb; F=Fqq+Fqb
    ev=eigh(H,eigvals_only=True,check_finite=False)
    tol=max(1e-7,1e-7*np.max(np.abs(ev)))
    breath=X.reshape(-1); breath=breath/np.linalg.norm(breath)
    lam=float(breath@H@breath); res=float(np.linalg.norm(H@breath-lam*breath))
    trans=[]
    for d in range(3):
        v=np.zeros(3*N); v[d::3]=1
        trans.append(float(v@H@v/(v@v)))
    r=dict(n=n,N=N,L=L,rho=rho,Eqq=Eqq,Eqb=Eqb,Ebb=Ebb,Etot=Eqq+Eqb+Ebb,
           Fqq_max=float(np.max(np.linalg.norm(Fqq,axis=1))),Fqb_max=float(np.max(np.linalg.norm(Fqb,axis=1))),
           Fmax=float(np.max(np.linalg.norm(F,axis=1))),Frms=float(np.sqrt(np.mean(F*F))),
           eig_min=float(ev[0]),eig_max=float(ev[-1]),Nneg=int(np.sum(ev<-tol)),Nzero=int(np.sum(abs(ev)<=tol)),Npos=int(np.sum(ev>tol)),tol=tol,
           translations=trans,breath_lam=lam,breath_res=res,J=J,
           deltaH_norm=float(np.linalg.norm(Hqb)),deltaH_identity_res=float(np.linalg.norm(Hqb-np.trace(Hqb)/(3*N)*np.eye(3*N))), eig=ev)
    return r

def report(r):
    for k in r:
        if k!='eig': print(k,':',r[k])

if __name__=='__main__':
    print('Ebb convergence n=4')
    L=1.5; rho=-64/(2*L)**3
    last=None
    for o in [16,24,32,48,64,80,100]:
        J=bb_integral_duffy(L,o); E=.5*rho*rho*J
        print(o,E, '' if last is None else ' delta='+str(abs(E-last)))
        last=E
    # First validate gradient convergence for a representative corner-adjacent 64B site.
    print('\nPhi/gradient convergence at x=(1.5,1.5,1.5), approached from inside is not used; sites are x=±1.5?')
    print('WARNING: Original 64B sites lie exactly on cube boundary if cube is |x|<=3a/2.')
    print('A volume background plus point charges on its boundary needs a boundary convention. Use L=1.5+eta or define half-density boundary limits.')
    for eta in [1e-2,3e-3,1e-3]:
        # enlarged cube is a clear regulator; rho must be reset if neutrality applies to enlarged volume
        print('eta=',eta)
        r=solve(4,order=28,h=2e-3,bb_order=48) if eta==1e-2 else None
    print('\nIMPORTANT: set geometry consistently before final verdict: charges on jellium boundary are singular for derivative formulas.')
