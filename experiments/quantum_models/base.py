import numpy as np, pandas as pd, json, itertools
from dataset_io import resolve_dataset_path
from scipy.spatial.distance import cdist
rng = np.random.default_rng(42)
DATA_PATH = resolve_dataset_path(__file__)
d = pd.read_csv(DATA_PATH)
kingdom = d.Kingdom.values
cod = [c.replace('U','T') for c in d.columns[1:]]
X = d.iloc[:,1:].values.astype(float); X = X/X.sum(1,keepdims=True)
# genetic code
bases='TCAG'; aas='FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
code={a+b+c:aas[i] for i,(a,b,c) in enumerate(itertools.product(bases,repeat=3))}
aa=np.array([code[c] for c in cod])
deg=np.array([(aa==a).sum() for a in aa],float)
gc=np.array([[ch in 'GC' for ch in c] for c in cod],float)  # 64x3
def obs(P):
    H=-(P*np.log(P+1e-12)).sum(1)
    GC=P@gc.mean(1); GC1,GC2,GC3=[P@gc[:,k] for k in range(3)]
    Dsyn=P@np.log(deg)
    # ENC (Wright) from frequencies, no small-sample correction; sense codons only, ignoring stops
    F={}
    for a in set(aa)-{'*'}:
        m=aa==a; s=P[:,m].sum(1)+1e-12; p=P[:,m]/s[:,None]; F[a]=((p**2).sum(1),m.sum())
    def avgF(k): 
        v=[F[a][0] for a in F if F[a][1]==k]; return np.mean(v,0) if v else None
    ENC=2.0
    for k in (2,3,4,6):
        f=avgF(k); nk=sum(1 for a in F if F[a][1]==k)
        ENC=ENC+ nk/np.clip(f,1e-9,None) if k!=3 else ENC
    # k=3 (Ile) handled as avg of F2 and F4 per Wright
    ENC=ENC+1/np.clip((avgF(2)+avgF(4))/2,1e-9,None)
    ENC=np.clip(ENC,20,61)
    return np.column_stack([P,H,GC,GC1,GC2,GC3,Dsyn,ENC])
Psi=obs(X); names=cod+['H','GC','GC1','GC2','GC3','Dsyn','ENC']
GC3=Psi[:,names.index('GC3')]
# split by organism
idx=rng.permutation(len(X)); ntr=int(.8*len(X)); tr,te=idx[:ntr],idx[ntr:]
def pairs(ids, metric='jensenshannon', k=10):
    """directed pairs i->j: j = nearest neighbour (within kingdom, same split) among k-NN having higher GC3"""
    out=[]
    for kg in np.unique(kingdom):
        ii=ids[kingdom[ids]==kg]
        if len(ii)<3: continue
        D=cdist(X[ii],X[ii],metric=metric); np.fill_diagonal(D,np.inf)
        for a in range(len(ii)):
            nn=np.argsort(D[a])[:k]
            up=[b for b in nn if GC3[ii[b]]>GC3[ii[a]]]
            if up: out.append((ii[a],ii[up[0]]))
    return np.array(out)
def fit(A,B,lam):
    return np.linalg.solve(A.T@A+lam*np.eye(A.shape[1]),A.T@B)
def zs(Z,mu,sd): return (Z-mu)/sd
def run(metric='jensenshannon',lam=1e-3,sub=None,label=''):
    Ptr=pairs(tr,metric); Pte=pairs(te,metric)
    mu=Psi[tr].mean(0); sd=Psi[tr].std(0)+1e-12
    Z=zs(Psi,mu,sd)
    ext=slice(0,Z.shape[1]); cd=slice(0,64)
    res={'n_train_pairs':len(Ptr),'n_test_pairs':len(Pte)}
    Y=X[Pte[:,1]]
    def rmse(Yh): return float(np.sqrt(((Y-Yh)**2).sum(1).mean()))
    res['RMSE_persistence']=rmse(X[Pte[:,0]])
    km={k:X[tr][kingdom[tr]==k].mean(0) for k in np.unique(kingdom)}
    res['RMSE_kingdom_mean']=rmse(np.array([km[kingdom[i]] for i in Pte[:,0]]))
    # codon-only linear (EDMD with dictionary = x)
    Kc=fit(Z[Ptr[:,0]][:,cd],Z[Ptr[:,1]][:,cd],lam)
    Yh=(Z[Pte[:,0]][:,cd]@Kc)*sd[cd]+mu[cd]; res['RMSE_koopman_codon']=rmse(Yh)
    Ke=fit(Z[Ptr[:,0]],Z[Ptr[:,1]],lam)
    Yh=(Z[Pte[:,0]]@Ke)[:,cd]*sd[cd]+mu[cd]; res['RMSE_koopman_extended']=rmse(Yh)
    # residual (drift) model: K regularised toward identity
    Ke2=fit(Z[Ptr[:,0]],Z[Ptr[:,1]]-Z[Ptr[:,0]],lam)+np.eye(Z.shape[1])
    Yh=(Z[Pte[:,0]]@Ke2)[:,cd]*sd[cd]+mu[cd]; res['RMSE_koopman_extended_toward_identity']=rmse(Yh)
    # paired bootstrap test extended vs codon (per-pair squared errors)
    e_c=((Y-((Z[Pte[:,0]][:,cd]@Kc)*sd[cd]+mu[cd]))**2).sum(1)
    e_e=((Y-((Z[Pte[:,0]]@Ke)[:,cd]*sd[cd]+mu[cd]))**2).sum(1)
    diff=e_c-e_e; bs=[rng.choice(diff,len(diff)).mean() for _ in range(2000)]
    res['H2_meandiff_MSE(codon-ext)']=float(diff.mean()); res['H2_CI95']=[float(np.percentile(bs,2.5)),float(np.percentile(bs,97.5))]
    return res,Ptr,Pte,Z,Ke,Ke2
def spec(K):
    w=np.linalg.eigvals(K); o=np.argsort(-np.abs(w)); return w[o]
