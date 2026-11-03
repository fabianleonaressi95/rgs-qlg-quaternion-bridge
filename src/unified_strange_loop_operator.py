"""Unified finite strange-loop operator: Cantor/Parry, Q6, Bach, Escher, Godel.

The operator is a finite combinatorial/spectral model. Its only certified
intertwiner to the 69,120-node PRISM Laplacian is the constant-mode map;
nontrivial spectral equivalence is deliberately left as an open problem.
"""
from __future__ import annotations
from hashlib import sha256
import json
import numpy as np

PSI=1.465571231876768
N_PRISM=69120


def bach_matrix():
    return np.array([[0.,0.,1.],[1.,0.,0.],[0.,1.,1.]])


def q6_laplacian():
    A=np.zeros((64,64))
    for x in range(64):
        for k in range(6):
            A[x,x^(1<<k)]=1.
    return 6*np.eye(64)-A


def gray2(i):
    i=int(i)%4
    return i^(i>>1)


def torus_q6_embedding():
    ids=[]
    for i in range(4):
        for j in range(4):
            gi,gj=gray2(i),gray2(j)
            bits=((gi>>1)&1,gi&1,(gj>>1)&1,gj&1,0,0)
            ids.append(sum(int(b)<<k for k,b in enumerate(bits)))
    if len(set(ids))!=16:
        raise ArithmeticError("Gray embedding is not injective")
    return np.asarray(ids,dtype=int)


def torus_laplacian():
    L=np.zeros((16,16))
    for i in range(4):
        for j in range(4):
            u=4*i+j
            for ni,nj in (((i+1)%4,j),((i-1)%4,j),(i,(j+1)%4),(i,(j-1)%4)):
                v=4*ni+nj
                L[u,u]+=1.
                L[u,v]-=1.
    return L


def escher_lifted_torus_laplacian():
    E=np.zeros((64,16))
    ids=torus_q6_embedding()
    E[ids,np.arange(16)]=1.
    return E@torus_laplacian()@E.T


def cantor_parry_q6_laplacian():
    """Symmetrized stationary-flow Laplacian coupling Cantor emissions to Q6."""
    A=bach_matrix()
    r=np.array([1.,PSI**-1,PSI])
    P=A*r[None,:]/(PSI*r[:,None])
    vals,vecs=np.linalg.eig(P.T)
    k=int(np.argmin(np.abs(vals-1.)))
    pi=np.real(vecs[:,k])
    if pi.sum()<0: pi=-pi
    pi=pi/pi.sum()
    phi=(1.+np.sqrt(5.))/2.
    qphi=np.exp(-1./phi**2)/phi
    q=qphi*r/float(np.prod(r)**(1./3.))
    p1=q/(1.+q); p0=1./(1.+q)
    F=np.zeros((192,192))
    for i in range(3):
        for x in range(64):
            u=64*i+x
            for j in range(3):
                if A[i,j]==0: continue
                bit_index=(i+j)%6
                for b,pb in ((0,p0[i]),(1,p1[i])):
                    y=x^(int(b)<<bit_index)
                    F[u,64*j+y]+=pi[i]*P[i,j]*pb/64.
    W=.5*(F+F.T)
    L=np.diag(W.sum(axis=1))-W
    return L,{"parry_transition":P,"stationary":pi,"cantor_p0":p0,"cantor_p1":p1}


def build_strange_loop_operator(kappa=1.0, sheet_weight=1.0):
    """Symmetric operator on Parry(3) x Cantor-prefix Q6(64) x sheets(2)."""
    if kappa<0 or sheet_weight<0:
        raise ValueError("weights must be nonnegative")
    Lcq,measure=cantor_parry_q6_laplacian()
    Lt=escher_lifted_torus_laplacian()
    Lsheet=np.array([[1.,-1.],[-1.,1.]])
    H=(np.kron(Lcq,np.eye(2))
       +kappa*np.kron(np.kron(np.eye(3),Lt),np.eye(2))
       +sheet_weight*np.kron(np.eye(192),Lsheet))
    return H,{"L_cantor_parry_q6":Lcq,"measure":measure,"L_torus_lifted":Lt,"L_sheet":Lsheet}



def doppler_moire_phase_audit(
    omega0=1.0, beta=0.1, k1=(1.0,0.0), k2=(1.0,0.08),
    projection=None, phase0=0.0, t=1.0
):
    """Doppler temporal generator and Moire spatial phase on the 64 Q6 vertices.

    projection maps six Q6 bit coordinates to a physical 2D plane (2x6).
    The default is an explicitly declared coordinate chart, not a derived
    RGS metric. Frequencies and wavevectors retain their input units.
    """
    if not (0.0 <= abs(beta) < 1.0):
        raise ValueError("beta must satisfy |beta|<1")
    k1=np.asarray(k1,dtype=float); k2=np.asarray(k2,dtype=float)
    if k1.shape!=(2,) or k2.shape!=(2,):
        raise ValueError("k1 and k2 must be two-dimensional wavevectors")
    if projection is None:
        projection=np.zeros((2,6)); projection[0,0]=1.; projection[1,1]=1.
    P=np.asarray(projection,dtype=float)
    if P.shape!=(2,6):
        raise ValueError("projection must have shape (2,6)")
    dk=k2-k1
    positions=np.asarray([P@np.array([(x>>a)&1 for a in range(6)],dtype=float)
                          for x in range(64)])
    phase=positions@dk+float(phase0)
    moire_potential=1.-np.cos(phase)
    # Observer moves along the first physical axis; c=1. Transform both
    # wave four-vectors, then form the *difference* of observed frequencies.
    gamma=1./np.sqrt(1.-beta*beta)
    omega1_obs=gamma*(float(omega0)-float(beta)*float(k1[0]))
    omega2_obs=gamma*(float(omega0)-float(beta)*float(k2[0]))
    delta_omega=omega2_obs-omega1_obs
    khat=k1/np.linalg.norm(k1) if np.linalg.norm(k1)>0 else np.array([1.,0.])
    doppler_factor=gamma*(1.-beta*float(khat[0]))
    # Keep the temporal generator separate from the dimensionless spatial
    # Laplacian; the relative phase is Delta omega*t - Delta k dot x.
    temporal_phase=delta_omega*float(t)
    spacetime_phase=temporal_phase-phase
    coherence=np.exp(-1j*spacetime_phase)
    return {
        "delta_k":dk.tolist(),
        "moire_length":float(2.*np.pi/np.linalg.norm(dk)) if np.linalg.norm(dk)>0 else None,
        "moire_phase":phase.tolist(),
        "moire_potential":moire_potential.tolist(),
        "doppler_factor":float(doppler_factor),
        "omega1_observed":float(omega1_obs),
        "omega2_observed":float(omega2_obs),
        "delta_omega_observed":float(delta_omega),
        "beat_period":float(2.*np.pi/abs(delta_omega)) if abs(delta_omega)>0 else None,
        "temporal_phase":float(temporal_phase),
        "spacetime_phase":spacetime_phase.tolist(),
        "coherence_real":coherence.real.tolist(),
        "coherence_imag":coherence.imag.tolist(),
        "phase_identity":"Delta theta(x,t)=Delta omega*t-Delta k dot x+phase0; Doppler is temporal, Moire is spatial",
        "projection":P.tolist(),
        "caveat":"projection and beta are model inputs; the RGS vierbein/action must determine them for a physical prediction",
    }


def build_doppler_moire_operator(kappa=1.0,sheet_weight=1.0,
                                  moire_strength=0.0,omega0=1.0,beta=0.0,
                                  k1=(1.0,0.0),k2=(1.0,0.08),projection=None):
    """Return static Hermitian spatial operator and separate Doppler generator."""
    if moire_strength<0:
        raise ValueError('moire_strength must be nonnegative')
    H,parts=build_strange_loop_operator(kappa,sheet_weight)
    audit=doppler_moire_phase_audit(omega0=omega0,beta=beta,k1=k1,k2=k2,projection=projection)
    vm=np.asarray(audit["moire_potential"])
    V=np.tile(vm,3)
    V=np.repeat(V,2)
    H_ext=H+float(moire_strength)*np.diag(V)
    omega=np.full(384,audit["delta_omega_observed"],dtype=float)
    return H_ext,omega,{**parts,"doppler_moire":audit,"moire_diagonal":V}


def _prism_constant_mode_intertwiner(H):
    """Factorized rank-one T=v u*: maps the constant mode to PRISM constant mode."""
    n=H.shape[0]
    u=np.ones(n)/np.sqrt(n)
    # PRISM has 69,120 vertices; normalized constant vector is represented by
    # its scalar normalization, avoiding allocation of a 69,120 x 384 matrix.
    v_scale=1./np.sqrt(N_PRISM)
    Hu=H@u
    source_res=float(np.linalg.norm(Hu))
    target_res=0.0  # L_PRISM 1=0 for every Cartesian product of cycles.
    return {
        "representation":"T = v u^T, factorized; v=1/sqrt(69120) * 1",
        "rank":1,
        "source_dimension":n,
        "target_dimension":N_PRISM,
        "source_kernel_residual":source_res,
        "target_kernel_residual":target_res,
        "intertwining_residual":float(source_res),
        "interpretation":"exact constant-mode intertwiner only; no nonconstant spectral equivalence",
    }


def audit_unified_operator(kappa=1.0,sheet_weight=1.0):
    H,parts=build_strange_loop_operator(kappa,sheet_weight)
    eig=np.linalg.eigvalsh(H)
    one=np.ones(H.shape[0])
    bach=bach_matrix()
    torus_ids=torus_q6_embedding()
    torus_edges=sum(1 for i in range(4) for j in range(4)
                    for ni,nj in (((i+1)%4,j),(i,(j+1)%4)))
    payload={
        "operator":"H=L_Cantor-Parry-Q6 (x) I_sheet + kappa I_Bach (x) L_Torus_Q6 (x) I_sheet + gamma I_192 (x) (I_sheet-J_sheet)",
        "shape":list(H.shape),"kappa":float(kappa),"sheet_weight":float(sheet_weight),
        "bach_matrix":bach.tolist(),"torus_embedding":torus_ids.tolist(),"cantor_branch":"golden state-dependent emissions on Parry transitions update Q6 bit (i+j) mod 6",
        "q6_vertices":64,"q6_edges":192,"torus_vertices":16,"torus_cycle_edges":torus_edges,
        "sheet_involution":"J(s)=1-s",
    }
    canonical=json.dumps(payload,sort_keys=True,separators=(",",":"))
    intertwiner=_prism_constant_mode_intertwiner(H)
    return {
        "audit":"UNIFIED_CANTOR_Q6_BACH_ESCHER_GODEL_STRANGE_LOOP",
        "status":"finite operator and rank-one PRISM constant-mode intertwiner verified; nontrivial intertwining open",
        "dimension":int(H.shape[0]),
        "symmetry_residual":float(np.linalg.norm(H-H.T)),
        "constant_mode_residual":float(np.linalg.norm(H@one)),
        "minimum_eigenvalue":float(eig[0]),
        "negative_eigenvalue_count":int(np.sum(eig < -1e-9)),
        "zero_mode_count":int(np.sum(np.abs(eig)<1e-9)),
        "bach_supergolden_polynomial_residual":float(PSI**3-PSI**2-1),
        "bach_matrix_identity_residual":float(np.linalg.norm(np.linalg.matrix_power(bach.astype(int),3)-np.linalg.matrix_power(bach.astype(int),2)-np.eye(3))),
        "torus_embedding_injective":len(set(torus_ids.tolist()))==16,
        "torus_lift_constant_residual":float(np.linalg.norm(parts["L_torus_lifted"]@np.ones(64))),
        "escher_sheet_involution_residual":float(np.linalg.norm(np.linalg.matrix_power(np.array([[0.,1.],[1.,0.]]),2)-np.eye(2))),
        "godel_finite_fingerprint_sha256":sha256(canonical.encode()).hexdigest(),
        "godel_caveat":"finite canonical fingerprint, not a Gödel incompleteness proof",
        "prism_intertwiner":intertwiner,
        "operator_definition":payload["operator"],
        "open_obligations":[
            "replace the rank-one constant-mode intertwiner with a nontrivial invariant-subspace map",
            "compare the 384-dimensional operator with the actual 69,120-node PRISM operator beyond its zero mode",
            "derive kappa and sheet_weight from the PSU action rather than treating them as model parameters",
            "derive the Bach substitution from a score-level symbolic coding",
            "construct a nontrivial Escher Z2 cocycle and its lifted cover",
            "prove convergence of finite Cantor cylinders to an infinite transfer operator",
            "define a formal arithmetic theory before making any Gödel incompleteness claim",
        ],
    }


def barycentric_gradient_moire_audit(
    positions, weights, delta_k, metric=None, delta_omega=0.0,
    phase0=0.0, t=0.0, eta=1.0
):
    """Exact barycentric/phase audit plus second-order Moire expansion.

    positions: (N,d) coordinate array; weights: nonnegative length-N weights.
    delta_k is a covector in these coordinates. metric is a constant SPD
    covariant metric g_ab; identity is used only when metric is omitted.
    """
    X=np.asarray(positions,dtype=float)
    w=np.asarray(weights,dtype=float)
    dk=np.asarray(delta_k,dtype=float)
    if X.ndim!=2 or w.shape!=(X.shape[0],) or dk.shape!=(X.shape[1],):
        raise ValueError("incompatible positions, weights, or delta_k shapes")
    if not np.all(np.isfinite(X)) or not np.all(np.isfinite(w)) or not np.all(np.isfinite(dk)):
        raise ValueError("inputs must be finite")
    if np.any(w<0) or w.sum()<=0:
        raise ValueError("weights must be nonnegative with positive sum")
    if eta<0:
        raise ValueError("eta must be nonnegative")
    w=w/w.sum()
    d=X.shape[1]
    g=np.eye(d) if metric is None else np.asarray(metric,dtype=float)
    if g.shape!=(d,d) or not np.allclose(g,g.T,atol=1e-12):
        raise ValueError("metric must be a symmetric (d,d) matrix")
    if np.linalg.eigvalsh(g)[0]<=0:
        raise ValueError("metric must be positive definite")
    xb=w@X
    Y=X-xb
    cov=(Y*w[:,None]).T@Y
    phase_b=float(delta_omega)*float(t)-float(dk@xb)+float(phase0)
    q=Y@dk
    phases=phase_b-q
    exact=float(np.sum(w*(1.-np.cos(phases))))
    char=np.sum(w*np.exp(-1j*q))
    exact_characteristic=float(1.-np.real(np.exp(1j*phase_b)*char))
    quadratic=float(1.-np.cos(phase_b)+.5*np.cos(phase_b)*(dk@cov@dk))
    dtheta=-dk
    grad_theta=np.linalg.solve(g,dtheta)
    # Force is minus the metric gradient of eta*(1-cos(theta)).
    grad_potential=float(eta)*np.sin(phases)[:,None]*grad_theta[None,:]
    force=-grad_potential
    return {
        "barycenter":xb.tolist(),
        "centered_positions":Y.tolist(),
        "weight_sum":float(w.sum()),
        "first_moment_residual":float(np.linalg.norm(w@Y)),
        "covariance":cov.tolist(),
        "covariance_eigenvalues":np.linalg.eigvalsh(cov).tolist(),
        "phase_at_barycenter":phase_b,
        "relative_phases":phases.tolist(),
        "characteristic_real":float(char.real),
        "characteristic_imag":float(char.imag),
        "moire_average_exact":exact,
        "moire_average_characteristic":exact_characteristic,
        "moire_average_quadratic":quadratic,
        "moire_quadratic_index":float(dk@cov@dk),
        "phase_covector_gradient":dtheta.tolist(),
        "metric_phase_gradient":grad_theta.tolist(),
        "moire_force":force.tolist(),
        "metric":g.tolist(),
        "expansion_order":"second order in delta_k dot centered_position; exact characteristic also returned",
        "caveat":"constant metric in the supplied coordinate chart; a varying RGS metric requires covariant transport and a specified connection",
    }


def moving_barycenter_audit(positions, weights, velocities, weight_rates=None):
    """Differentiate X_B=sum(w_i X_i), including changing normalized weights."""
    X=np.asarray(positions,dtype=float)
    w=np.asarray(weights,dtype=float)
    V=np.asarray(velocities,dtype=float)
    if X.ndim!=2 or V.shape!=X.shape or w.shape!=(X.shape[0],):
        raise ValueError("incompatible positions, velocities, or weights")
    if np.any(w<0) or w.sum()<=0:
        raise ValueError("weights must be nonnegative with positive sum")
    w=w/w.sum()
    wd=np.zeros_like(w) if weight_rates is None else np.asarray(weight_rates,dtype=float)
    if wd.shape!=w.shape:
        raise ValueError("weight_rates must have shape (N,)")
    # Rates are required to be tangent to the normalized simplex.
    wd=wd-wd.sum()*w
    xb=w@X
    y=X-xb
    xbdot=w@V+wd@X
    identity_residual=float(np.linalg.norm(xbdot-(w@V+wd@y)))
    return {
        "barycenter":xb.tolist(),
        "barycenter_velocity":xbdot.tolist(),
        "normalized_weight_rate_sum":float(wd.sum()),
        "centered_weight_rate_term":(wd@y).tolist(),
        "derivative_identity_residual":identity_residual,
    }
