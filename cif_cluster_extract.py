#!/usr/bin/env python3
"""CIF -> cluster candidates -> 2D tiling vertices for Al-Pd-Mn Xi audit.
Requires: pip install numpy pymatgen scipy
"""
import argparse, csv, hashlib, json, math, sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
try:
    from pymatgen.core import Structure
    from scipy.spatial import ConvexHull
except ImportError as e:
    raise SystemExit('Dipendenza mancante. Esegui: pip install numpy pymatgen scipy\n'+str(e))

TAU=(1+math.sqrt(5))/2

def sha256(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for q in iter(lambda:f.read(1<<20),b''): h.update(q)
    return h.hexdigest()

def species(site):
    return str(site.specie) if site.is_ordered else str(site.species)

def wrapped_delta(frac_a, frac_b):
    d=np.asarray(frac_b)-np.asarray(frac_a)
    return d-np.rint(d)

def cart_delta(struct, i, j):
    return wrapped_delta(struct[i].frac_coords,struct[j].frac_coords) @ struct.lattice.matrix

def pca_plane(points):
    c=points.mean(axis=0)
    _,_,vh=np.linalg.svd(points-c,full_matrices=False)
    u,v=vh[0],vh[1]
    return c,u,v

def project_plane(points,c,u,v):
    d=points-c
    return np.column_stack((d@u,d@v))

def local_neighborhood(struct, center_idx, rmax):
    center=struct[center_idx]
    neighbors=struct.get_neighbors(center,rmax)
    data=[]
    for nb in neighbors:
        try: j=int(nb.index)
        except Exception: continue
        if j==center_idx: continue
        vec=cart_delta(struct,center_idx,j)
        data.append((j,float(np.linalg.norm(vec)),vec,species(struct[j])))
    return sorted(data,key=lambda q:q[1])

def cluster_score(struct, i, r_core, r_shell):
    ns=local_neighborhood(struct,i,r_shell)
    core=[q for q in ns if q[1]<=r_core]
    shell=[q for q in ns if r_core<q[1]<=r_shell]
    # Heuristica trasparente: coordinazione del core vicina a 12 e popolazione shell significativa.
    score=max(0.0,1.0-abs(len(core)-12)/12)*0.55 + min(len(shell)/30,1.0)*0.30
    if species(struct[i]).startswith('Pd') or species(struct[i]).startswith('Mn'): score+=0.15
    return min(score,1.0),core,shell

def find_clusters(struct,r_core,r_shell,threshold):
    out=[]
    for i in range(len(struct)):
        score,core,shell=cluster_score(struct,i,r_core,r_shell)
        if score>=threshold:
            out.append({'center_index':i,'species':species(struct[i]),'frac_coords':[float(x) for x in struct[i].frac_coords],
                        'cart_coords_A':[float(x) for x in struct[i].coords],'score':float(score),
                        'core_neighbors':len(core),'shell_neighbors':len(shell),
                        'core_species':dict(Counter(q[3] for q in core)),'shell_species':dict(Counter(q[3] for q in shell))})
    return sorted(out,key=lambda x:(-x['score'],x['center_index']))

def select_layer(clusters,axis,thickness):
    pts=np.array([x['cart_coords_A'] for x in clusters])
    t=pts@axis
    median=float(np.median(t))
    keep=[c for c,val in zip(clusters,t) if abs(val-median)<=thickness/2]
    return keep,median

def make_vertices(clusters,axis):
    pts=np.array([x['cart_coords_A'] for x in clusters])
    c,u,v=pca_plane(pts)
    xy=project_plane(pts,c,u,v)/10.0
    return xy,{'origin_A':c.tolist(),'u':u.tolist(),'v':v.tolist()}

def main():
    ap=argparse.ArgumentParser(description='Estrae candidati cluster e vertici 2D da CIF.')
    ap.add_argument('--cif',required=True)
    ap.add_argument('--outdir',default='cif_cluster_output')
    ap.add_argument('--r-core',type=float,default=3.5,help='raggio core in Angstrom')
    ap.add_argument('--r-shell',type=float,default=8.0,help='raggio shell in Angstrom')
    ap.add_argument('--threshold',type=float,default=.55)
    ap.add_argument('--axis',choices=['a','b','c'],default='b',help='asse delle colonne presunte')
    ap.add_argument('--layer-thickness',type=float,default=2.0,help='spessore slice in Angstrom')
    args=ap.parse_args()
    if args.r_shell<=args.r_core: ap.error('--r-shell deve essere > --r-core')
    struct=Structure.from_file(args.cif)
    out=Path(args.outdir); out.mkdir(parents=True,exist_ok=True)
    clusters=find_clusters(struct,args.r_core,args.r_shell,args.threshold)
    if len(clusters)<3: raise SystemExit('Meno di tre candidati: amplia r-shell, riduci threshold o controlla il CIF.')
    axis={'a':struct.lattice.matrix[0],'b':struct.lattice.matrix[1],'c':struct.lattice.matrix[2]}[args.axis]
    axis=np.asarray(axis)/np.linalg.norm(axis)
    layer,layer_center=select_layer(clusters,axis,args.layer_thickness)
    if len(layer)<3:
        layer=clusters
        layer_center=None
    xy,plane=make_vertices(layer,axis)
    with open(out/'cluster_candidates.csv','w',newline='') as f:
        fields=['center_index','species','score','core_neighbors','shell_neighbors','frac_coords','cart_coords_A','core_species','shell_species']
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader()
        for x in clusters: w.writerow({k:json.dumps(x[k]) if isinstance(x[k],(dict,list)) else x[k] for k in fields})
    with open(out/'tiling_vertices.csv','w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['x_nm','y_nm','center_index','species','score']); w.writeheader()
        for q,x in zip(xy,layer): w.writerow({'x_nm':q[0],'y_nm':q[1],'center_index':x['center_index'],'species':x['species'],'score':x['score']})
    report={'engine':'CIF Cluster Extractor','schema_version':'1.0.0','created_utc':datetime.now(timezone.utc).isoformat(),
      'input':{'cif':str(args.cif),'sha256':sha256(args.cif),'sites':len(struct),'formula':struct.composition.reduced_formula,
               'lattice_A':struct.lattice.matrix.tolist()},
      'parameters':{'r_core_A':args.r_core,'r_shell_A':args.r_shell,'threshold':args.threshold,'column_axis':args.axis,'layer_thickness_A':args.layer_thickness},
      'result':{'cluster_candidates':len(clusters),'selected_layer':len(layer),'layer_center_projection_A':layer_center,
                'vertices_csv':str(out/'tiling_vertices.csv'),'clusters_csv':str(out/'cluster_candidates.csv'),'plane':plane},
      'warning':'I candidati sono identificati con una euristica di coordinazione, non costituiscono un riconoscimento automatico definitivo di cluster Mackay.'}
    (out/'cluster_extraction_report.json').write_text(json.dumps(report,indent=2,ensure_ascii=False))
    print(f'CIF: {args.cif}\nFormula: {struct.composition.reduced_formula}; siti: {len(struct)}\nCandidati: {len(clusters)}; layer: {len(layer)}\nVertici: {out / "tiling_vertices.csv"}\nReport: {out / "cluster_extraction_report.json"}')
if __name__=='__main__': main()
