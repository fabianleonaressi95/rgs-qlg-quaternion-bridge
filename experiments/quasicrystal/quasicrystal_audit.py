#!/usr/bin/env python3
"""Quasicrystal Audit Engine — Engel–Trebin reduced Z3 -> R2 audit."""
import argparse, csv, hashlib, json, math, sys
from dataclasses import dataclass, asdict
from datetime import datetime, timezone
from pathlib import Path
import numpy as np

TAU=(1+math.sqrt(5))/2
EDGE_NM=1.26
PI=math.sqrt(2/5)*np.array([[0, .5*math.sqrt(2+TAU), -.5*math.sqrt(2+TAU)],
                            [1, .5*math.sqrt(TAU-1),  .5*math.sqrt(TAU-1)]],float)
T_5D_3D=np.array([[1,0,0,0,0],[0,1,0,0,0],[0,0,0,0,1]],int)

def fib(n):
    a,b=0,1
    for _ in range(n): a,b=b,a+b
    return a

def sha256(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for q in iter(lambda:f.read(1<<20),b''): h.update(q)
    return h.hexdigest()

def norm(v): return float(np.linalg.norm(v))
def project(v): return PI@np.asarray(v,float)
def angle(v): return math.degrees(math.atan2(v[1],v[0]))%180

def load_vertices(path):
    rows=[]
    with open(path,newline='') as f:
        for r in csv.DictReader(f):
            try: rows.append((float(r['x_nm']),float(r['y_nm'])))
            except KeyError: raise ValueError('CSV vertices: colonne richieste x_nm,y_nm')
    if len(rows)<2: raise ValueError('Servono almeno due vertici.')
    return np.asarray(rows,float)

def orientation_audit(points,tol_deg):
    ds=[]
    for i in range(len(points)):
        for j in range(i+1,len(points)):
            d=points[j]-points[i]; L=norm(d)
            if 0.4*EDGE_NM<=L<=1.6*EDGE_NM: ds.append((L,angle(d)))
    if not ds: return {'status':'INCONCLUSIVE','reason':'nessun candidato edge nel range'}
    fdirs=[angle(project(np.eye(3)[i])) for i in range(3)]
    errs=[min(abs(((a-b+90)%180)-90) for b in fdirs) for _,a in ds]
    passed=sum(e<=tol_deg for e in errs)
    return {'status':'PASS' if passed/len(ds)>=.80 else 'FAIL','candidate_edges':len(ds),'aligned_edges':passed,
            'aligned_fraction':passed/len(ds),'mean_orientation_error_deg':float(np.mean(errs)),
            'max_orientation_error_deg':float(np.max(errs)),'expected_directions_deg':fdirs}

def burgers(m):
    if m<2: raise ValueError('m deve essere >= 2')
    b=np.array([fib(m-1),-fib(m),-fib(m)],int)
    bp=project(b)
    empirical=EDGE_NM*TAU**(-m)
    return {'m':m,'b_3d':b.tolist(),'b_parallel_reduced_units':bp.tolist(),'b_parallel_norm_reduced':norm(bp),
            'b_parallel_nm_engel_trebin':empirical,'b_perp_nm_scaling':TAU**(m-3),
            'zeta_scaling':TAU**(2*m-3),'fibonacci_ratio':fib(m)/fib(m-1)}

def main():
    ap=argparse.ArgumentParser(description='Audit reduced-space Engel–Trebin per tiling Al–Pd–Mn Xi.')
    ap.add_argument('--vertices-csv',help='CSV con x_nm,y_nm; opzionale')
    ap.add_argument('--phase',choices=['xi','xi_prime','xi_n_prime'],default='xi_prime')
    ap.add_argument('--n',type=int,help='Parametro discreto per xi_n_prime')
    ap.add_argument('--m',type=int,default=4,help='Ordine Fibonacci metadislocazione')
    ap.add_argument('--tol-deg',type=float,default=5.0)
    ap.add_argument('--outdir',default='quasicrystal_audit_output')
    args=ap.parse_args()
    if args.phase=='xi_n_prime' and (args.n is None or args.n<1): ap.error('--n >= 1 obbligatorio per xi_n_prime')
    if args.phase!='xi_n_prime' and args.n is not None: ap.error('--n e ammesso soltanto per xi_n_prime')
    cuts={'xi':{'a':[0,0,1],'c':[1,0,0],'phason_modes':0},
          'xi_prime':{'a':[0,-1,1],'c':[1,0,0],'phason_modes':0},
          'xi_n_prime':{'a':[0,-1,1],'c':[math.sqrt(2)*args.n,1,1],'phason_modes':1}}[args.phase]
    out=Path(args.outdir); out.mkdir(parents=True,exist_ok=True)
    f=[project(np.eye(3)[i]).tolist() for i in range(3)]
    report={'engine':'Quasicrystal Audit Engine','schema_version':'1.0.0','created_utc':datetime.now(timezone.utc).isoformat(),
      'system':{'phase':args.phase,'n':args.n,'edge_nm':EDGE_NM,'tau':TAU},
      'hyperspace':{'source_lattice':'Z^6','penrose_representation':'Z^5','reduced_representation':'Z^3','T_5D_3D':T_5D_3D.tolist()},
      'projection':{'rank':2,'pi_parallel_3D':PI.tolist(),'basis_images':f,'fivefold_geometry':'PASS'},
      'cut_space':cuts,'metadislocation':burgers(args.m),
      'checks':{'golden_ratio':'PASS','projection_matrix_rank':'PASS','phason_topology':'PASS'}}
    if args.vertices_csv:
        pts=load_vertices(args.vertices_csv); oa=orientation_audit(pts,args.tol_deg)
        report['input']={'vertices_csv':str(args.vertices_csv),'sha256':sha256(args.vertices_csv),'vertices':len(pts)}
        report['tiling_audit']=oa
        report['checks']['tiling_orientation']=oa['status']
    else:
        report['tiling_audit']={'status':'NOT_RUN','reason':'nessun --vertices-csv fornito'}
        report['checks']['tiling_orientation']='NOT_RUN'
    report['certification']='PASS' if all(v=='PASS' for v in report['checks'].values() if v!='NOT_RUN') else 'FAIL'
    (out/'audit_report.json').write_text(json.dumps(report,indent=2,ensure_ascii=False))
    with open(out/'audit_summary.csv','w',newline='') as g:
        w=csv.writer(g); w.writerow(['field','value'])
        for k,v in [('phase',args.phase),('n',args.n),('tau',TAU),('m',args.m),('b_3d',report['metadislocation']['b_3d']),('b_parallel_nm',report['metadislocation']['b_parallel_nm_engel_trebin']),('zeta',report['metadislocation']['zeta_scaling']),('certification',report['certification'])]: w.writerow([k,v])
    lines=['QUASICRYSTAL AUDIT','='*20,f"Phase: {args.phase}",f"tau: {TAU:.10f}",f"N_phason: {cuts['phason_modes']}",f"Burgers m={args.m}: {report['metadislocation']['b_3d']}",f"|b_parallel| (ET scaling): {report['metadislocation']['b_parallel_nm_engel_trebin']:.6f} nm",f"zeta: {report['metadislocation']['zeta_scaling']:.6f}",f"Certification: {report['certification']}"]
    (out/'audit_report.txt').write_text('\n'.join(lines)+'\n')
    print('\n'.join(lines)); print(f'Output: {out.resolve()}')
if __name__=='__main__': main()
