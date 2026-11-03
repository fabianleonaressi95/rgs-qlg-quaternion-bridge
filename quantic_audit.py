#!/usr/bin/env python3
"""Quantic Module — local momentum audit from electron density grid."""
import argparse, json, math, sys
from datetime import datetime, timezone
from pathlib import Path
import numpy as np

def load_density_grid(path):
    """Carica un CSV con x,y,z,rho oppure un formato semplice."""
    xs,ys,zs,rhos=[],[],[],[]
    with open(path,newline='') as f:
        import csv
        for r in csv.DictReader(f):
            xs.append(float(r['x'])); ys.append(float(r['y'])); zs.append(float(r['z'])); rhos.append(float(r['rho']))
    if not xs: raise ValueError('Nessun dato caricato')
    xs=np.asarray(xs); ys=np.asarray(ys); zs=np.asarray(zs); rhos=np.asarray(rhos)
    # Assumiamo griglia regolare; deduciamo dimensioni
    nx=len(np.unique(xs)); ny=len(np.unique(ys)); nz=len(np.unique(zs))
    rho=rhos.reshape(nz,ny,nx)
    x=np.unique(xs); y=np.unique(ys); z=np.unique(zs)
    dx=x[1]-x[0] if len(x)>1 else 1.0
    dy=y[1]-y[0] if len(y)>1 else 1.0
    dz=z[1]-z[0] if len(z)>1 else 1.0
    return {'x':x,'y':y,'z':z,'rho':rho,'dx':dx,'dy':dy,'dz':dz}

def grad_rho(rho,dx,dy,dz):
    gy,gx,gz=np.gradient(rho,dx,dy,dz)
    return gx,gy,gz

def local_momentum_bohm(rho,gx,gy,gz):
    """p_loc = -ħ/2 * ∇ρ/ρ (Bohmian, forma semplificata)."""
    hbar=1.054571817e-34
    rho_safe=np.where(rho>1e-12,rho,1e-12)
    px=-hbar/2*gx/rho_safe
    py=-hbar/2*gy/rho_safe
    pz=-hbar/2*gz/rho_safe
    return px,py,pz

def main():
    ap=argparse.ArgumentParser(description='Audit quantistico: momento locale da densità elettronica.')
    ap.add_argument('--density-csv',required=True,help='CSV con x,y,z,rho (griglia regolare)')
    ap.add_argument('--outdir',default='quantic_audit_output')
    args=ap.parse_args()
    out=Path(args.outdir); out.mkdir(parents=True,exist_ok=True)
    data=load_density_grid(args.density_csv)
    gx,gy,gz=grad_rho(data['rho'],data['dx'],data['dy'],data['dz'])
    px,py,pz=local_momentum_bohm(data['rho'],gx,gy,gz)
    p_norm=np.sqrt(px**2+py**2+pz**2)
    stats={'mean_p':float(np.mean(p_norm)),'max_p':float(np.max(p_norm)),'min_p':float(np.min(p_norm)),
           'std_p':float(np.std(p_norm)),'grid_points':int(p_norm.size)}
    report={'engine':'Quantic Module','schema_version':'1.0.0','created_utc':datetime.now(timezone.utc).isoformat(),
            'input':{'density_csv':str(args.density_csv),'grid_shape':list(data['rho'].shape)},
            'method':'Bohmian local momentum p_loc = -ħ/2 ∇ρ/ρ',
            'statistics':stats}
    (out/'quantic_report.json').write_text(json.dumps(report,indent=2,ensure_ascii=False))
    lines=['QUANTIC AUDIT (Local Momentum)','='*35,f"Input: {args.density_csv}",f"Grid: {data['rho'].shape}",
           f"Mean |p_loc|: {stats['mean_p']:.6e} kg·m/s",f"Max |p_loc|: {stats['max_p']:.6e} kg·m/s",
           f"Std |p_loc|: {stats['std_p']:.6e} kg·m/s"]
    (out/'quantic_report.txt').write_text('\n'.join(lines)+'\n')
    print('\n'.join(lines)); print(f'Output: {out.resolve()}')
if __name__=='__main__': main()
