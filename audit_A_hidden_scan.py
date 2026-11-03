"""Hidden-sector requirement vs beta4 (uses trajectory_A.csv; F=T; M_f^2=1, m4 Mg2=1). Finite differences => ~1e-6 accuracy."""
import numpy as np
d=np.genfromtxt("trajectory_A.csv",delimiter=",",names=True)
t,H,Y,X,Hm,r=d['t'],d['H'],d['Y'],d['X'],d['H'],d['rho']; Hm=d['H']
Hm=d['H_']if 'H_' in d.dtype.names else Hm
Hminus=H/Y
print("prompt's condition H+ - X H- : max|.| =",np.max(np.abs(H-X*Hminus)),"(should vanish if (b1+2b2Y+b3Y^2)!=0; it does not)")
print("correct condition H+ - Y H- : max|.| =",np.max(np.abs(H-Y*Hminus)))
m=np.abs(t)<6
for b4 in (30.0,10.0,3.0,0.0):
    Vm=2-3*Y+b4*Y**3
    rho=3*Hminus**2-Vm/Y**3
    drho=np.gradient(rho,t)
    p=-rho-drho/(3*Hminus*X)*np.where(np.abs(Hminus)<1e-9,np.nan,1)  # N-=X, 3 N- H- = 3 a-'/a-
    # use rho_dot = -3 (a-'/a-)(rho+p),  a-'/a- = X H-
    ok=np.abs(Hminus)>1e-3
    print("beta4=%5.1f: rho_- in [%.3f, %.3f] ; (rho+p)_- in [%.3f, %.3f] (|H-|>1e-3, |t|<6)"%(b4,rho[m].min(),rho[m].max(),(rho+p)[m&ok].min(),(rho+p)[m&ok].max()))
