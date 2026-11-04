import numpy as np, json
exec(open('base.py').read())
k=10; lam=1e-3
mu=Psi[tr].mean(0); sd=Psi[tr].std(0)+1e-12; Z=(Psi-mu)/sd; cd=slice(0,64)
def knn_graph(ids):
    G={}
    for kg in np.unique(kingdom):
        ii=ids[kingdom[ids]==kg]
        if len(ii)<3: continue
        D=cdist(X[ii],X[ii],metric='jensenshannon'); np.fill_diagonal(D,np.inf)
        for a in range(len(ii)): G[ii[a]]=ii[np.argsort(D[a])[:k]]
    return G
Gtr=knn_graph(tr); Gte=knn_graph(te)
def directed(G):   # real rule: nearest neighbour with higher GC3
    return np.array([(a,[b for b in nb if GC3[b]>GC3[a]][0]) for a,nb in G.items() if any(GC3[b]>GC3[a] for b in nb)])
Ptr=directed(Gtr); Pte=directed(Gte)
Y=X[Pte[:,1]]
def stats(P):
    K=fit(Z[P[:,0]],Z[P[:,1]],lam); w=np.linalg.eigvals(K); w=w[np.argsort(-np.abs(w))]
    Yh=(Z[Pte[:,0]]@K)[:,cd]*sd[cd]+mu[cd]
    return dict(m1=abs(w[0]),m2=abs(w[1]),m3=abs(w[2]),m5=abs(w[4]),
        maxim=np.abs(w[:10].imag).max(), ncx=int((np.abs(w[:10].imag)>1e-9).sum()),
        asym=np.linalg.norm(K-K.T)/np.linalg.norm(K), rmse=float(np.sqrt(((Y-Yh)**2).sum(1).mean())),
        drift=float((Z[P[:,1]]-Z[P[:,0]])[:,names.index('GC3')].mean()))
real=stats(Ptr)
B=300
def nullrun(gen):
    return [stats(gen()) for _ in range(B)]
def flip():
    P=Ptr.copy(); m=rng.random(len(P))<.5; P[m]=P[m][:,::-1]; return P
def randnb():
    return np.array([(a,rng.choice(nb)) for a,nb in Gtr.items()])
def randnb_matched():  # same sources as real pairs, random neighbour from kNN (graph-preserving, direction-free)
    return np.array([(a,rng.choice(Gtr[a])) for a in Ptr[:,0]])
def gc3_shuffled():  # direction rule applied to permuted GC3 labels: keeps graph and rule form, breaks biology
    perm=rng.permutation(len(GC3)); g=GC3[perm]
    return np.array([(a,[b for b in nb if g[b]>g[a]][0]) for a,nb in Gtr.items() if any(g[b]>g[a] for b in nb)])
nulls={'flip_direction':nullrun(flip),'random_neighbour':nullrun(randnb_matched),'gc3_label_permuted':nullrun(gc3_shuffled)}
out={'real':real,'n_pairs':len(Ptr)}
for n,L in nulls.items():
    o={}
    for key in real:
        v=np.array([x[key] for x in L]); 
        o[key]={'null_mean':float(v.mean()),'null_sd':float(v.std()),'z':float((real[key]-v.mean())/(v.std()+1e-12)),
                'p_two_sided':float((1+(np.abs(v-v.mean())>=abs(real[key]-v.mean())).sum())/(B+1))}
    out[n]=o
# bootstrap real for sd
bs=[stats(Ptr[rng.integers(0,len(Ptr),len(Ptr))]) for _ in range(200)]
out['real_boot_sd']={key:float(np.std([b[key] for b in bs])) for key in real}
json.dump(out,open('results_graphnull.json','w'),indent=1)
print('REAL',{a:round(float(b),5) for a,b in real.items()})
print('bootsd',{a:round(b,5) for a,b in out['real_boot_sd'].items()})
for n in nulls:
    print('\n',n)
    for key in real: print(f"  {key:7s} null {out[n][key]['null_mean']:.5f}±{out[n][key]['null_sd']:.5f}  z={out[n][key]['z']:.1f}  p={out[n][key]['p_two_sided']:.3f}")
