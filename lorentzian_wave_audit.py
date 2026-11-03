"""Lorentzian convergence prototype for a 1+1 wave equation.
This is not a GR solver. It tests the separate hyperbolic mechanism.
Equation: u_tt = c^2 u_xx - gamma u_t + F(t,x).
Periodic domain, leapfrog/Newmark-style centered scheme.
"""
import numpy as np
from dataclasses import dataclass

@dataclass
class Result:
    dx: float
    dt: float
    l2_error: float
    energy: np.ndarray
    constraint: np.ndarray

def exact(x, t, c=1.0):
    return np.cos(2*np.pi*(x-c*t))

def run(nx, T=1.0, c=1.0, cfl=0.4, gamma=0.0):
    L=1.0; dx=L/nx; dt=cfl*dx/c; nt=int(np.ceil(T/dt)); dt=T/nt
    x=np.arange(nx)*dx
    u=exact(x,0,c); v=2*np.pi*c*np.sin(2*np.pi*x)
    lap=lambda z:(np.roll(z,-1)-2*z+np.roll(z,1))/dx**2
    # Taylor start, second order in time
    u_prev=u.copy(); u=u+dt*v+0.5*dt**2*(c*c*lap(u)-gamma*v)
    energies=[]; constraints=[]
    for n in range(1,nt):
        vhalf=(u-u_prev)/dt
        u_next=2*u-u_prev+dt**2*(c*c*lap(u)-gamma*vhalf)
        vnext=(u_next-u_prev)/(2*dt)
        grad=(np.roll(u,-1)-np.roll(u,1))/(2*dx)
        E=0.5*np.mean(vnext**2+c*c*grad**2)
        # discrete PDE residual at current time
        res=(u_next-2*u+u_prev)/dt**2-c*c*lap(u)+gamma*vhalf
        energies.append(E); constraints.append(np.sqrt(np.mean(res**2)))
        u_prev,u=u,u_next
    err=np.sqrt(np.mean((u-exact(x,T,c))**2))
    return Result(dx,dt,err,np.array(energies),np.array(constraints))

def convergence_table(levels=(50,100,200,400)):
    rows=[]
    for nx in levels:
        r=run(nx)
        rows.append((r.dx,r.dt,r.l2_error,r.energy[-1],r.constraint.max()))
    for i in range(len(rows)-1):
        p=np.log(rows[i][2]/rows[i+1][2])/np.log(rows[i][0]/rows[i+1][0])
        rows[i]+=(p,)
    return rows

if __name__ == "__main__":
    for row in convergence_table(): print(row)
