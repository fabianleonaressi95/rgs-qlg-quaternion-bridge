import numpy as np, pandas as pd, json, itertools
from scipy.spatial.distance import cdist
rng = np.random.default_rng(42)
d = pd.read_csv('codon_train.csv')
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
res,Ptr,Pte,Z,Ke,Ke2=run()
w=spec(Ke); w2=spec(Ke2)
res['top_eig_modulus_K']=[float(x) for x in np.abs(w[:8])]
res['top_eig_arg_K']=[float(x) for x in np.angle(w[:8])]
res['n_complex_in_top10']=int((np.abs(w[:10].imag)>1e-9).sum())
res['top_eig_modulus_K_identityreg']=[float(x) for x in np.abs(w2[:8])]
# bootstrap over training pairs
B=200; boots=[]
for b in range(B):
    s=rng.integers(0,len(Ptr),len(Ptr)); K=fit(Z[Ptr[s,0]],Z[Ptr[s,1]],1e-3); boots.append(np.abs(spec(K)[:8]))
boots=np.array(boots)
res['boot_mean_mod']=boots.mean(0).tolist(); res['boot_sd_mod']=boots.std(0).tolist()
# null: permute targets within kingdom (destroys pair relation, keeps marginals & group sizes)
nul=[]
for b in range(B):
    T=Ptr[:,1].copy()
    for kg in np.unique(kingdom):
        m=kingdom[Ptr[:,0]]==kg
        T[m]=rng.permutation(T[m])
    K=fit(Z[Ptr[:,0]],Z[T],1e-3); nul.append(np.abs(spec(K)[:8]))
nul=np.array(nul)
res['null_mean_mod']=nul.mean(0).tolist(); res['null_max_mod']=nul.max(0).tolist()
res['H1_lambda2_real']=float(np.abs(w[1])); res['H1_lambda2_null_95pct']=float(np.percentile(nul[:,1],95))
res['H1_p']=float((1+(nul[:,1]>=np.abs(w[1])).sum())/(B+1))
# H4: proximity of eigenvalue moduli to 1/phi and 2 (|lam|<1 expected)
phi_inv=(np.sqrt(5)-1)/2
def mind(v): return np.min(np.abs(v-phi_inv))
real_d=mind(np.abs(w)); null_d=[]
for b in range(B):
    T=Ptr[:,1].copy()
    for kg in np.unique(kingdom):
        m=kingdom[Ptr[:,0]]==kg; T[m]=rng.permutation(T[m])
    null_d.append(mind(np.abs(spec(fit(Z[Ptr[:,0]],Z[T],1e-3)))))
res['H4_min_|mod-1/phi|_real']=float(real_d); res['H4_null_median']=float(np.median(null_d)); res['H4_p_closer_by_chance']=float((1+(np.array(null_d)<=real_d).sum())/(B+1))
res['max_modulus_overall']=float(np.abs(w).max())
res['tau_top(pair steps)']=[float(-1/np.log(min(r,0.999999))) for r in np.abs(w[:5])]
# robustness of pairing: alternative distance / lambda
rob={}
for metric,lam in [('euclidean',1e-3),('jensenshannon',1e-1),('jensenshannon',1e-5)]:
    r,_,_,_,K_,_=run(metric,lam); ww=spec(K_)
    rob[f'{metric}_lam{lam}']={'RMSE_persist':r['RMSE_persistence'],'RMSE_codon':r['RMSE_koopman_codon'],'RMSE_ext':r['RMSE_koopman_extended'],'RMSE_ext_toward_id':r['RMSE_koopman_extended_toward_identity'],'top_mod':[float(x) for x in np.abs(ww[:5])]}
res['robustness']=rob
# static geometry control: per-kingdom summary
res['gc3_by_kingdom']={k:[float(GC3[kingdom==k].mean()),float(GC3[kingdom==k].std())] for k in np.unique(kingdom)}
res['ENC_by_kingdom']={k:float(Psi[kingdom==k,-1].mean()) for k in np.unique(kingdom)}
json.dump(res,open('results.json','w'),indent=1)
print(json.dumps(res,indent=1))
