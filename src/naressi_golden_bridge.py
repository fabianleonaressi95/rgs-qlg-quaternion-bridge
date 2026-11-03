"""Naressi golden-coupling bridge: QAT, TFIM and Q6."""
from __future__ import annotations
from math import sqrt
import numpy as np
from scipy.linalg import expm
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.linalg import eigsh

PHI = (1.0 + sqrt(5.0)) / 2.0

def golden_parameters(phi=PHI, alpha=1.0, n=1):
    phi, alpha = float(phi), float(alpha)
    if not np.isfinite(phi) or phi <= 1.0:
        raise ValueError("phi must be finite and greater than one")
    if not np.isfinite(alpha) or alpha <= 0.0:
        raise ValueError("alpha must be finite and positive")
    if not isinstance(n, (int, np.integer)) or n < 1:
        raise ValueError("n must be a positive integer")
    theta = 2.0*np.pi/(phi*phi)
    q = np.exp(-1.0/(phi*phi))/phi
    j = alpha*q**int(n)
    return {"phi":phi,"alpha":alpha,"n":float(n),"q_phi":float(q),
            "theta_phi":float(theta),"J_n":float(j),
            "chi_n":float(j*np.cos(theta)),
            "omega_photon":float(phi*phi),"omega_phonon":float(phi)}

def golden_decay_to_q6_tau(phi=PHI):
    if not np.isfinite(phi) or phi <= 1.0:
        raise ValueError("phi must be finite and greater than one")
    q = np.exp(-1.0/(phi*phi))/phi
    return float(-np.log(q))

def q6_entropy(tau):
    tau=float(tau)
    if not np.isfinite(tau) or tau < 0:
        raise ValueError("tau must be finite and non-negative")
    p=1.0/(1.0+np.exp(2.0*tau))
    if p == 0.0:
        return 0.0
    return float(-6.0*(p*np.log(p)+(1.0-p)*np.log1p(-p)))

def golden_q6_bridge(phi=PHI):
    tau=golden_decay_to_q6_tau(phi)
    q=float(np.exp(-tau))
    return {"dictionary":"q_phi = exp(-tau), tau = -log(q_phi)",
            "status":"explicit spectral ansatz; geometric derivation remains open",
            "phi":float(phi),"q_phi":q,"tau":tau,"q6_entropy_nats":q6_entropy(tau)}

def _annihilation(d):
    if not isinstance(d,(int,np.integer)) or d<2:
        raise ValueError("Fock dimension must be an integer >= 2")
    a=np.zeros((d,d),dtype=complex)
    for n in range(1,d): a[n-1,n]=np.sqrt(n)
    return a

def qat_hamiltonian(dimension=3,phi=PHI,alpha=1.0,n_site=1):
    """Hermitian finite-Fock model: wp*n_a+wf*n_b+chi*(n_a*n_b+b+b†)."""
    p=golden_parameters(phi,alpha,n_site)
    a0=_annihilation(dimension); ident=np.eye(dimension,dtype=complex)
    a=np.kron(a0,ident); b=np.kron(ident,a0)
    na=a.conj().T@a; nb=b.conj().T@b
    h=(p["omega_photon"]*na+p["omega_phonon"]*nb
       +p["chi_n"]*(na@nb+b+b.conj().T))
    return 0.5*(h+h.conj().T)

def _pure_bipartite_entropy(state,d):
    s=np.linalg.svd(np.asarray(state).reshape(d,d),compute_uv=False)
    p=s*s; p=p[p>1e-14]
    return float(-np.sum(p*np.log(p)))

def audit_qat_unitarity(dimension=3,times=(0.,.1,.25,.5,1.),
                        phi=PHI,alpha=1.,n_site=1):
    h=qat_hamiltonian(dimension,phi,alpha,n_site)
    herm=float(np.linalg.norm(h-h.conj().T))
    vals,vecs=np.linalg.eigh(h)
    psi0=np.zeros(dimension*dimension,dtype=complex)
    psi0[0]=1/np.sqrt(2); psi0[dimension]=1/np.sqrt(2)
    rows=[]; max_u=0.; max_n=0.
    for t in times:
        u=(vecs*np.exp(-1j*vals*float(t)))@vecs.conj().T
        ur=float(np.linalg.norm(u.conj().T@u-np.eye(dimension**2)))
        psi=u@psi0; nr=abs(float(np.vdot(psi,psi).real)-1.)
        rows.append({"time":float(t),"entropy_nats":_pure_bipartite_entropy(psi,dimension),
                     "unitarity_residual":ur,"norm_residual":nr})
        max_u=max(max_u,ur); max_n=max(max_n,nr)
    return {"audit":"QAT_GOLDEN_TWO_MODE_UNITARITY","dimension_per_mode":dimension,
            "hermiticity_residual":herm,"max_unitarity_residual":max_u,
            "max_state_norm_residual":max_n,"rows":rows,
            "pass":bool(herm<1e-12 and max_u<1e-10 and max_n<1e-10)}

def tfim_golden_hamiltonian(n_sites,g=1.,alpha=1.,phi=PHI):
    """Open TFIM CSR matrix, J_n=alpha*q_phi**n, big-endian site indexing."""
    if not isinstance(n_sites,(int,np.integer)) or n_sites<2:
        raise ValueError("n_sites must be an integer >= 2")
    if n_sites>24: raise ValueError("CSR exact diagonalization is limited to N <= 24")
    if not np.isfinite(g): raise ValueError("g must be finite")
    p=golden_parameters(phi,alpha,1); dim=1<<n_sites
    states=np.arange(dim,dtype=np.int64)
    diag=np.zeros(dim,dtype=float)
    for site in range(n_sites):
        bit=(states>>(n_sites-1-site))&1
        diag-=float(g)*(1.-2.*bit)
    rows=[states]; cols=[states]; data=[diag]
    for site in range(n_sites-1):
        mask=(1<<(n_sites-1-site))|(1<<(n_sites-2-site))
        rows.append(states); cols.append(states^mask)
        data.append(np.full(dim,-float(alpha)*p["q_phi"]**(site+1)))
    return coo_matrix((np.concatenate(data),(np.concatenate(rows),np.concatenate(cols))),
                      shape=(dim,dim),dtype=float).tocsr()

def tfim_ground_state_audit(n_sites,g=1.,alpha=1.,phi=PHI):
    h=tfim_golden_hamiltonian(n_sites,g,alpha,phi)
    vals,vecs=eigsh(h,k=2,which="SA",tol=1e-11)
    order=np.argsort(vals); vals=vals[order]; vec=vecs[:,order[0]]
    psi=vec.reshape(1<<(n_sites//2),1<<(n_sites-n_sites//2))
    s=np.linalg.svd(psi,compute_uv=False); p=s*s; p=p[p>1e-14]
    return {"n_sites":int(n_sites),"ground_energy":float(vals[0]),
            "first_excited_energy":float(vals[1]),"gap":float(vals[1]-vals[0]),
            "half_chain_entropy_nats":float(-np.sum(p*np.log(p))),
            "hermiticity_residual":float(np.linalg.norm((h-h.T).data))}

def audit_golden_naressi_bridge(phi=PHI):
    pars=golden_parameters(phi,1.,1); bridge=golden_q6_bridge(phi)
    return {"audit":"NARESSI_GOLDEN_QAT_TFIM_Q6_BRIDGE","parameters":pars,
            "q6_spectral_dictionary":bridge,
            "evidence":{"golden_decay_ratio":"exact from declared coupling law",
              "chi_projection":"exact scalar definition chi=J_1*cos(theta_phi)",
              "q6_entropy":"exact finite-Q6 formula evaluated at tau=-log(q_phi)",
              "tau_to_geometry":"OPEN: derive from QAT coframe/action",
              "tfim_to_q6_state_map":"OPEN: TFIM ground state is not automatically Q6 heat-kernel state"}}
