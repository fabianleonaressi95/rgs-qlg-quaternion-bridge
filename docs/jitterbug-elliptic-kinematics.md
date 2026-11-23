# Exact symmetric Jitterbug kinematics

This reconstruction distinguishes the vertex ellipse parameter theta from
the rigid triangular-face rotation gamma.

## Vertex formula

For each cuboctahedron vertex with zero coordinate k and nonzero signs
s_i,s_j, set

u=(s_i e_i+s_j e_j)/sqrt(2), v=e_k cross u, sigma=-s_i s_j.

Then the vertex trajectory is

p(theta)=Z(theta)u+sigma Y(theta)v,

Z=cos(theta)/sqrt(3-2cos(theta)^2),
Y=sin(theta)/sqrt(3-2cos(theta)^2).

This gives 3Y^2+Z^2=1. At theta=0 the cuboctahedron is recovered;
theta=+/-pi/4 gives the octahedral endpoints.

## Quaternionic rigid-face pose

For a triangular face indexed by epsilon=(epsilon_1,epsilon_2,epsilon_3),
define u_epsilon=epsilon/sqrt(3), and

gamma(theta)=atan(sqrt(3) tan(theta)),
s_epsilon=-epsilon_1 epsilon_2 epsilon_3.

The scalar-first Hamilton quaternion is

q_epsilon(theta)=(cos(gamma/2),
 s_epsilon u_epsilon sin(gamma/2)).

The face translation is

t_epsilon(theta)=sqrt(2/3)(cos(gamma)-1)u_epsilon.

The rigid pose acts actively on column vectors:

X(theta)=R(q_epsilon)X(0)+t_epsilon.

## Numerical audit

The script checks all 24 edge lengths, all eight rigid triangles, and shared
joint coincidence over 101 parameter samples on [-pi/4,pi/4].
Run: python scripts/jitterbug_elliptic_kinematics.py

The formula is consistent with the ellipse/rotation-angle derivation in
Robert W. Gray, The Jitterbug Motion (2002):
https://rwgrayprojects.com/rbfnotes/polyhed/jitterbug/paper01/jbCurve.pdf

This is a symmetric one-parameter path through the bar-and-joint framework;
it is not a claim that the unconstrained point-hinge assembly has only one
degree of freedom.