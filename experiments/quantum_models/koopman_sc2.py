import numpy as np, pandas as pd, json
rng=np.random.default_rng(7)
df=pd.read_csv('sc2_codon_counts.csv',index_col=0)
meta=pd.read_csv('ncov/data/example_metadata.tsv',sep='\t').set_index('strain')
ok=df.ncodons>=9000; df=df[ok]
codons=list(df.columns[:-1]); C=df[codons].values; X=C/C.sum(1,keepdims=True)
date=pd.to_datetime(meta.loc[df.index,'date'])
week=((date-pd.Timestamp('2019-12-23')).dt.days//7).values   # weeks from Mon 2019-12-23
# keep weeks with >=8 genomes
wk,cn=np.unique(week,return_counts=True); keep=wk[cn>=8]
print('genomes',len(df),'weeks kept',dict(zip(wk[np.isin(wk,keep)],cn[np.isin(wk,keep)])))
mask=np.isin(week,keep); X=X[mask]; week=week[mask]; W=np.array(sorted(set(week)))
def weekly(X,week,idx=None):
    return np.array([X[week==w].mean(0) for w in W])
Xw=weekly(X,week)
n_w=np.array([(week==w).sum() for w in W])
# 1) signal vs sampling noise
dobs=np.mean(np.sum(np.diff(Xw,axis=0)**2,1))
noise=[]
for t in range(len(W)-1):
    v=lambda w: X[week==w].var(0,ddof=1).sum()/ (week==w).sum()
    noise.append(v(W[t])+v(W[t+1]))
dnoise=np.mean(noise)
print('mean squared week-to-week change observed %.3e, expected from sampling noise alone %.3e'%(dobs,dnoise))
# 2) Koopman (EDMD, ridge) on PCA-reduced state of weekly means
def pca_state(Xw,r):
    mu=Xw.mean(0); U,S,Vt=np.linalg.svd(Xw-mu,full_matrices=False)
    return (Xw-mu)@Vt[:r].T, mu, Vt[:r], S
def edmd(Z,lam):
    A,B=Z[:-1],Z[1:]; return np.linalg.solve(A.T@A+lam*np.eye(A.shape[1]),A.T@B)
def forward_chain(Z,lam,start=6):
    pe=[];pp=[];pm=[]
    for t in range(start,len(Z)-1):
        K=edmd(Z[:t+1],lam); pe.append(Z[t]@K-Z[t+1]); pp.append(Z[t]-Z[t+1]); pm.append(Z[:t+1].mean(0)-Z[t+1])
    f=lambda a: float(np.sqrt(np.mean(np.sum(np.array(a)**2,1))))
    return f(pe),f(pp),f(pm)
_,_,_,S=pca_state(Xw,3); print('PCA singular values',np.round(S[:6],5))
out={'n_genomes':int(len(df)),'weeks':[int(w) for w in W],'n_per_week':n_w.tolist(),'msd_week_change_obs':float(dobs),'msd_week_change_noise':float(dnoise)}
res={}
for r in (1,2,3):
    Z,mu,V,_=pca_state(Xw,r); lam=1e-8*np.trace(Z.T@Z)
    K=edmd(Z,lam); w=np.linalg.eigvals(K); w=w[np.argsort(-np.abs(w))]
    e,p,m=forward_chain(Z,lam)
    # bootstrap genomes within week
    bm=[];
    for b in range(300):
        Xb=np.array([X[week==wv][rng.integers(0,(week==wv).sum(),(week==wv).sum())].mean(0) for wv in W])
        Zb=(Xb-mu)@V.T; Kb=edmd(Zb,lam); wb=np.linalg.eigvals(Kb); bm.append(np.sort(np.abs(wb))[::-1])
    bm=np.array(bm)
    # null: shuffle the week order of the weekly means
    nm=[];nr=[]
    for b in range(500):
        pi=rng.permutation(len(Z)); Zs=Z[pi]; Ks=edmd(Zs,lam); wn=np.linalg.eigvals(Ks); nm.append(np.sort(np.abs(wn))[::-1])
        nr.append(forward_chain(Zs,lam)[0])
    nm=np.array(nm); nr=np.array(nr)
    res[r]={'eig_real':[[float(z.real),float(z.imag)] for z in w],'modulus':np.abs(w).tolist(),
      'boot_mod_mean':bm.mean(0).tolist(),'boot_mod_sd':bm.std(0).tolist(),
      'null_mod_mean':nm.mean(0).tolist(),'null_mod_p95':np.percentile(nm,95,axis=0).tolist(),
      'p_lambda1_ge_real':float((1+(nm[:,0]>=abs(w[0])).sum())/501),
      'RMSE_koopman':e,'RMSE_persistence':p,'RMSE_mean':m,
      'RMSE_koopman_null_mean':float(nr.mean()),'RMSE_koopman_null_p5':float(np.percentile(nr,5)),
      'p_rmse_le_real':float((1+(nr<=e).sum())/501)}
    print(f"\nr={r}: |lambda|={np.round(np.abs(w),3)}  eig={np.round(w,3)}")
    print(f"  boot sd |l1|={bm[:,0].std():.3f}; null mean |l1|={nm[:,0].mean():.3f} p95={np.percentile(nm[:,0],95):.3f}  p={res[r]['p_lambda1_ge_real']:.3f}")
    print(f"  fwd-chain RMSE koopman={e:.2e} persistence={p:.2e} mean={m:.2e}  null-koopman mean={nr.mean():.2e}  p(rmse<=real)={res[r]['p_rmse_le_real']:.3f}")
out['koopman']=res
# 3) monotone drift check: does PC1 trend with time?
Z1=pca_state(Xw,1)[0][:,0]; rho=np.corrcoef(W,Z1)[0,1]; print('\nPC1 vs week corr',round(rho,3))
out['pc1_time_corr']=float(rho)
# top drifting codons early (wk<=8) vs late
early=X[week<=W[2]].mean(0); late=X[week>=W[-2]].mean(0)
d=(late-early); o=np.argsort(-np.abs(d))[:6]; print('largest codon shifts', [(codons[i],float(d[i])) for i in o])
out['top_shifts']=[(codons[i],float(d[i])) for i in o]
json.dump(out,open('results_sc2.json','w'),indent=1)
