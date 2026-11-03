import numpy as np, json, hashlib
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.linalg import eigsh, LinearOperator
from math import comb
mods=np.array([64,6,6,6,5],dtype=np.int64); D=5; V=int(np.prod(mods)); E=D*V; P=10; F2=P*V
strides=np.array([np.prod(mods[d+1:]) for d in range(D)],dtype=np.int64)
verts=np.arange(V,dtype=np.int64)
# all shifted targets: shape D,V
shifts=np.empty((D,V),dtype=np.int64)
for d in range(D):
 s=strides[d]; m=mods[d]; q=(verts//s)%m
 shifts[d]=verts + np.where(q<m-1,s,-(m-1)*s)
# B1 vectorized
src=np.tile(verts,D); dst=shifts.ravel(order='C'); edges=np.repeat(np.arange(E,dtype=np.int64),2)
rows=np.empty(2*E,dtype=np.int32); cols=np.empty(2*E,dtype=np.int32); dat=np.empty(2*E,dtype=np.int8)
rows[0::2]=src; rows[1::2]=dst; cols[0::2]=np.repeat(np.arange(E),1); cols[1::2]=np.repeat(np.arange(E),1); dat[0::2]=-1; dat[1::2]=1
B1=coo_matrix((dat,(rows,cols)),shape=(V,E)).tocsr()
# B2
rs=[]; cs=[]; ds=[]; fbase=0
for i in range(D):
 for j in range(i+1,D):
  vi=shifts[i]; vj=shifts[j];
  ei_v=i*V+verts; ej_vi=j*V+vi; ei_vj=i*V+vj; ej_v=j*V+verts
  fs=np.arange(fbase,fbase+V,dtype=np.int64)
  rs.append(np.concatenate([ei_v,ej_vi,ei_vj,ej_v])); cs.append(np.concatenate([fs,fs,fs,fs])); ds.append(np.concatenate([np.ones(V,dtype=np.int8),np.ones(V,dtype=np.int8),-np.ones(V,dtype=np.int8),-np.ones(V,dtype=np.int8)]))
  fbase+=V
B2=coo_matrix((np.concatenate(ds),(np.concatenate(rs),np.concatenate(cs))),shape=(E,F2)).tocsr()
B1B2=B1@B2
# exact scalar spectrum
arr=np.array([0.0])
for n in mods:
 a=np.array([2-2*np.cos(2*np.pi*k/n) for k in range(int(n))])
 arr=(arr[:,None]+a[None,:]).ravel()
arr.sort()
L0=B1@B1.T
first=arr[np.argmax(arr>1e-12)]
sha=hashlib.sha256()
for M in (B1,B2):
 sha.update(M.indptr.tobytes()); sha.update(M.indices.tobytes()); sha.update(M.data.tobytes())
report={'mods':mods.tolist(),'V':V,'E1':E,'E2':F2,'B1_shape':list(B1.shape),'B1_nnz':int(B1.nnz),'B2_shape':list(B2.shape),'B2_nnz':int(B2.nnz),'B1B2_nnz':int(B1B2.nnz),'graph_cycle_rank':E-V+1,'cells':[V,V*5,V*10,V*10,V*5,V],'betti':[1,5,10,10,5,1],'euler':0,'L0_exact_first12':arr[:12].tolist(),'L0_zero_mult':int(np.sum(arr<1e-12)),'L0_gap':float(first),'L0_gap_mult':int(np.sum(np.isclose(arr,first,atol=1e-12))),'L1_zero_mult_expected':5,'L1_gap_expected':float(first),'L1_first_nonzero_mult_expected':5*int(np.sum(np.isclose(arr,first,atol=1e-12))),'structural_sha256':sha.hexdigest(),'prism_recompute':False}
open('/mnt/data/hodge_69120_audit.json','w').write(json.dumps(report,indent=2))
np.savez_compressed('/mnt/data/hodge_69120_incidence.npz',B1_indptr=B1.indptr,B1_indices=B1.indices,B1_data=B1.data,B2_indptr=B2.indptr,B2_indices=B2.indices,B2_data=B2.data)
print(json.dumps(report,indent=2))
