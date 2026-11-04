#!/usr/bin/env python3
"""Genera CSV pilota di vertici per testare l'audit geometrico."""
import csv, math
from pathlib import Path
TAU=(1+math.sqrt(5))/2
EDGE=1.26
theta0=0
pts=[]
for k in range(5):
    theta=theta0+k*72*math.pi/180
    pts.append({'x_nm':EDGE*math.cos(theta),'y_nm':EDGE*math.sin(theta)})
    theta2=theta+36*math.pi/180
    pts.append({'x_nm':EDGE*TAU*math.cos(theta2),'y_nm':EDGE*TAU*math.sin(theta2)})
Path('output/test_vertices.csv').write_text('x_nm,y_nm\n'+''.join(f"{p['x_nm']:.6f},{p['y_nm']:.6f}\n" for p in pts))
print('output/test_vertices.csv')
