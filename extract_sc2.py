import gzip, numpy as np, pandas as pd, mappy, itertools, json
from Bio import SeqIO
ref=SeqIO.read('ncov/defaults/reference_seq.gb','genbank')
refseq=str(ref.seq).upper()
cds=[(f.qualifiers.get('gene',[f.qualifiers.get('product',['?'])[0]])[0],f) for f in ref.features if f.type=='CDS']
print([(n,len(f.location.extract(ref.seq))) for n,f in cds])
al=mappy.Aligner(seq=refseq,preset='asm5')
seqs={};k=None
for l in gzip.open('ncov/data/example_sequences.fasta.gz','rt'):
    if l[0]=='>': k=l[1:].split()[0]; seqs[k]=[]
    else: seqs[k].append(l.strip())
seqs={a:''.join(b).upper().replace('U','T') for a,b in seqs.items()}
def project(q):
    """return array of length len(ref): base aligned to each ref position ('-' if none)"""
    out=np.array(['-']*len(refseq),dtype='<U1')
    best=None
    for h in al.map(q):
        if best is None or h.mlen>best.mlen: best=h
    hits=sorted(al.map(q),key=lambda h:h.r_st)
    for h in hits:
        if h.mlen<500: continue
        qi=h.q_st; ri=h.r_st
        qs=q if h.strand==1 else mappy.revcomp(q)
        for n,op in h.cigar:
            if op in (0,7,8):
                out[ri:ri+n]=list(qs[qi:qi+n]); qi+=n; ri+=n
            elif op==1: qi+=n
            elif op in (2,3): ri+=n
            elif op in (4,5): pass
    return out
bases='TCAG'; codons=[''.join(c) for c in itertools.product(bases,repeat=3)]
cidx={c:i for i,c in enumerate(codons)}
rows=[];names=[]
for s,q in seqs.items():
    p=project(q)
    cnt=np.zeros(64); tot=0
    for n,f in cds:
        idx=np.array(f.location.extract(ref.seq) and [i for part in f.location.parts for i in (range(part.start,part.end) if part.strand!=-1 else range(part.end-1,part.start-1,-1))])
        b=''.join(p[idx])
        for i in range(0,len(b)-2,3):
            c=b[i:i+3]
            if c in cidx: cnt[cidx[c]]+=1
    rows.append(cnt); names.append(s)
C=np.array(rows); ntot=C.sum(1)
df=pd.DataFrame(C,columns=codons,index=names); df['ncodons']=ntot
df.to_csv('sc2_codon_counts.csv')
print(df.ncodons.describe())
