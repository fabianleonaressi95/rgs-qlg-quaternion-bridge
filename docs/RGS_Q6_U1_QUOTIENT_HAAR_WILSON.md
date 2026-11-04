# Q6 U(1) quotient Haar measure and Wilson action

Status: exact finite-lattice construction with numerical audits. It is not a
derivation of the fine-structure constant.

## 1. Add an explicit electromagnetic gauge candidate

Existing Q6 transport uses internal SU(2); that does not imply electromagnetism.
This contribution defines a separate compact U(1) link variable
(U_e=e^{i\theta_e}) on each of the 192 canonically oriented edges of Q6.
There are 64 vertices and 240 square faces. With incidence matrices
(D_0:C^0\to C^1), (D_1:C^1\to C^2), the exact identity is
(D_1D_0=0). Vertex phases act by
[
\theta' = \theta-D_0\lambda\pmod{2\pi},
\qquad F=D_1\theta\pmod{2\pi}.
]
Hence plaquette holonomies are gauge invariant.

## 2. Quotient and canonical finite measure

A deterministic spanning tree has 63 edges; the remaining 129 edges are
chords. Every gauge orbit has a unique tree gauge representative with tree
links equal to one, modulo the global U(1), which acts trivially. Thus
[
\mathcal A_{Q6}/\mathcal G_{Q6}\simeq U(1)^{192}/U(1)^{63}
\simeq U(1)^{129}.
]
The normalized quotient Haar probability measure is
[
d\mu_{Q6}=\prod_{e\in E\setminus T}\frac{d\theta_e}{2\pi},
\quad -\pi\leq\theta_e<\pi.
]
It is intrinsic to this compact abelian quotient and has total mass one.
Tree gauge fixing is a translation with unit Jacobian for Haar measure.

The restricted plaquette matrix (D_{1,C}) is (240\times129) and has rank
129. Since the full cubical 6-cell is contractible, there are no residual
flat link modes after quotienting by vertex gauge transformations. The
rank and singular values are audited in code.

## 3. Wilson action and Maxwell-type expansion

For bare (\beta>0),
[
S_{Q6}^{U(1)}[\theta]=\beta\sum_{f=1}^{240}
\left(1-\cos((D_1\theta)_f)\right).
]
It is exactly gauge invariant. At small flux,
[
S_{Q6}^{U(1)}=\frac{\beta}{2}\|D_1\theta\|^2
-\frac{\beta}{24}\sum_f(D_1\theta)_f^4+O(\|D_1\theta\|^6).
]
The leading term is a finite-lattice Maxwell-type quadratic form. The bare
coefficient beta is still free; it is not computed from a volume ratio.
The audit also verifies invariance under all 60 proper icosahedral actions.

## 4. Why this is not yet Wyler D5

The candidate (D_5\simeq SO_0(5,2)/(SO(5)\times SO(2))) has real dimension
10. This Q6 U(1) quotient has real dimension 129. It is therefore not D5,
and this construction supplies no dimension-preserving identification.
Its normalized compact Haar measure is not the invariant Bergman measure
of a noncompact bounded symmetric domain, whose total volume is infinite.

A future reduction must independently construct a 10-real-dimensional
effective sector/submanifold, its induced metric, and the Jacobian from the
Q6 measure. Selecting it to reproduce alpha would be fitting, not prediction.

## 5. Physical matching still open

Not supplied here: physical edge length; continuum gauge-potential and
electric-charge normalization; a four-dimensional continuum limit from the
finite six-cube; renormalization trajectory and matching to Thomson alpha;
or a derivation of the U(1) field from the RGS coframe/PSU variables.
Ultimately one needs
[
\Gamma[A]\supset\frac{1}{4g^2(\mu)}
\int d^4x\sqrt{|g|}F_{\mu\nu}F^{\mu\nu},
\qquad \alpha(\mu)=\frac{g^2(\mu)}{4\pi},
]
with every normalization fixed independently of alpha.

## Evidence ledger

**Exact by construction:** compact U(1) gauge action; 63+129 tree/chord
quotient; normalized product Haar measure; gauge-invariant Wilson action;
quadratic small-field expansion.

**Numerically audited:** (D_1D_0=0), reduced rank 129, gauge-orbit
representative, icosahedral invariance, small-field convergence.

**Open:** D5 sector and induced measure; scale/Jacobian; physical charge
normalization; continuum QED limit; prediction of alpha.

Reproduce:
```bash
PYTHONPATH=. python -m pytest tests/test_rgs_q6_u1_quotient.py -q
PYTHONPATH=. python -c "from src.rgs_q6_u1_quotient import audit_q6_u1_quotient; print(audit_q6_u1_quotient())"
```
