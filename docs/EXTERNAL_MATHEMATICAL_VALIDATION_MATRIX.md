# External Mathematical Validation Matrix

**Purpose:** turn external correspondence into falsifiable review questions. This register is not evidence of endorsement by any researcher.

Evidence labels: **PROVED** (proof under stated hypotheses); **NUMERICALLY VERIFIED** (reproducible finite computation); **CANDIDATE** (constructed, not established); **OPEN** (proof or independent validation missing).

## 1. Priority matrix

| ID | Object / claim | Internal basis | Independent review sought | Acceptance gate | Status |
|---|---|---|---|---|---|
| V1 | Q6 Wilson marginal and five-dimensional quotient | Integer incidence, saturated-kernel and fibre-sum constructions; Fourier–Haar marginal | Lattice gauge theory, integer lattices, harmonic analysis | Independently verify kernel saturation and rank-five quotient; derive the exact fibre theta sum; compare with a Gaussian only through Poisson summation plus a uniform remainder/tail bound | **OPEN** for uniform asymptotics and global measure push-forward. Rank five does not imply a type-IV domain or Bergman measure |
| V2 | Curved Ellis–RGS Mosco convergence | Curved G1 geometry and G4 common-L2 embedding; separate flat periodic P1 G6 theorem | FEEC/DEC, elliptic PDE, variational convergence | Uniform shape regularity and metric ellipticity; discrete Rellich compactness; liminf and recovery sequence in one fixed Hilbert space | **OPEN** for curved G5/G6. Flat-torus results cannot be transferred automatically |
| V3 | Dirac/Laplace spectra, kernel and gap | Finite 188-dimensional tangent audit and discrete operators | Spectral geometry, Dirac operators, heat kernels, operator convergence | Specify Hilbert spaces and self-adjoint domains; prove consistency/stability or resolvent convergence; control spectral projections, multiplicities, zero modes and gap | **NUMERICALLY VERIFIED** only at finite dimension; continuum interpretation **OPEN** |
| V4 | Quaternionic, hyperkähler or quaternion-Kähler interpretation | Finite quaternionic and tetrad constructions | Quaternionic geometry, spin geometry, twistor theory | Exhibit compatible global structures; verify integrability/parallelism and the relevant curvature conditions | **CANDIDATE**; algebraic quaternionic identities alone are insufficient |
| V5 | G2/Spin(7), calibrations and instantons | Discrete metric/gauge constructions | Exceptional holonomy and gauge theory | Construct defining differential forms; compute intrinsic torsion; prove torsion-free conditions or solve the instanton equation with analytic estimates | **OPEN**; finite symmetry sectors do not establish exceptional holonomy |
| V6 | Unified RGS–PSU / teleparallel / gauge action | Finite vierbein and teleparallel constructions | Geometric analysis, mathematical relativity, gauge theory | State one action and variables; derive Euler–Lagrange equations, symmetries, Bianchi identities and controlled limits | **OPEN**; no unification claim without this derivation |

## 2. Existing results: scope and limitations

- **Q6 finite algebra:** exact incidence identities and integer-rank calculations are finite combinatorial results, not a continuum limit.
- **Wilson fibre:** exact finite-lattice quotient and fibre sum do not justify replacing an integer theta sum by a Gaussian without a uniform error estimate.
- **Bergman/Hua:** a local Maxwell-chart pullback on its stated preimage is distinct from a global Wilson push-forward to Bergman measure.
- **Mosco:** the flat periodic P1 result applies to its regular flat mesh family. Curved Ellis–RGS G5/G6 needs its own compactness, liminf and recovery proof.
- **188D spectrum:** residuals, multiplicities and gaps describe the finite matrix under its conventions; they are not eigenvalue convergence or a continuum Dirac index theorem.
- **Quaternionic / exceptional geometry:** representation identities do not replace integrability, parallel differential forms, curvature identities or torsion-free conditions.

## 3. Routing external expertise

| Contacts | Primary review target | Do not ask them to infer |
|---|---|---|
| Ralf Hiptmair, Ragnar Winther, Douglas Arnold | V2: curved FEEC/DEC stability, compactness and Mosco | Endorsement of the whole RGS programme |
| Kaushik Kalyanaraman, Luis Wirth | V2/V3: independent DEC/FEEC implementation and benchmark interoperability | Continuum physics from code similarity |
| Luigi Del Debbio, Peter Boyle, Alessio Martini, Anita Maria Tabacco | V1: Wilson fibre, theta/Gaussian control, harmonic analysis and measure push-forward | Bergman measure from quotient dimension alone |
| Pierre Albin, Thomas Ourmières-Bonafos, Giovanni Landi | V3: Dirac operators, self-adjointness, spectral projections, index and convergence | A continuum mass gap from a finite matrix gap |
| Andrei Moroianu, Simon Salamon, Nobuhiro Honda | V4: spinors, twistor and quaternionic integrability | Hyperkähler geometry from quaternionic labels alone |
| Andrew J. Clarke, Spiro Karigiannis, Jason D. Lotay | V5: G2/Spin(7), torsion, calibrations and instantons | Exceptional holonomy from finite symmetry alone |
| Jeff Viaclovsky, Ioana Șuvaina, Jacopo Stoppa | V4/V6: elliptic geometric analysis, four-manifold constraints and Kähler compatibility | Existence from formal analogy |
| Alessandro Lucantonio | Applied benchmark in continuum mechanics and structure-preserving numerics | Number-theoretic or Bergman claims outside his expertise |

This is a suggested allocation of questions, not a claim that anyone has accepted, replied or endorsed the work. Verify current affiliations and contact details before follow-ups.

## 4. Correspondence-to-proof protocol

For each actual reply:
1. Record date, exact question, response, references and attribution/permission status.
2. Convert it into a precise lemma, counterexample, hypothesis list or benchmark—not a general endorsement.
3. Create one minimal repository artifact (proof note, test, benchmark or issue) and record its exact commit.
4. Seek an independent implementation or proof for computational claims.
5. Change an evidence label only when the relevant acceptance gate above is met.

Do not send another overlapping wave before consolidating replies. Technical priority: **V1 uniform lattice-sum error**, then **V2 curved Mosco G5/G6**, then **V3 operator and spectral convergence**.

## 5. Response tracker

No response is inferred from an email being drafted or sent. Update only from actual correspondence.

| Contact group | Topic | Status |
|---|---|---|
| Hiptmair, Winther, Arnold | Curved FEEC/DEC and Mosco | Not yet recorded |
| Del Debbio, Boyle, Martini, Tabacco | Wilson fibre and measure | Not yet recorded |
| Wirth, Kalyanaraman | DEC/FEEC implementation | Not yet recorded |
| Albin, Ourmières-Bonafos, Landi | Dirac and spectral convergence | Not yet recorded |
| Honda, Șuvaina, Stoppa, Viaclovsky, Moroianu, Cherkis, Clarke, Karigiannis, Lotay, Salamon | Special, spin and gauge geometry | Not yet recorded |
| Lucantonio | Computational mechanics benchmark | Not yet recorded |

**"Sent" is not "received", "reviewed" or "validated".**
