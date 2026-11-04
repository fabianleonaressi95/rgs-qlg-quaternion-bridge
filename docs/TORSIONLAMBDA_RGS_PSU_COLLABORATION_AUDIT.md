# Integration audit: Danny Lee Eldridge's TorsionLambda toolkit

**Purpose:** define a technically defensible collaboration interface between RGS–PSU–QLG and Danny Lee Eldridge's public TorsionLambda toolkit. This is an independent repository review, not a claim that either model has been observationally validated.

External project: https://github.com/dannyeldridge54/torsion-lambda-toolkit

Files inspected: README.md, docs/TECHNICAL_DERIVATION.md, src/models/torsion-lambda.ts, src/physics/consistency-checker.ts, docs/FALSIFICATION_CRITERIA.md, LICENSE.

## 1. Why this collaboration is relevant

TorsionLambda is positioned as a pre-simulation consistency and screening toolkit for a phenomenological Einstein–Cartan-inspired cosmological background. It provides a parameterized FLRW extension of LambdaCDM, background observables and limiting-case checks, energy-condition screening, and reproducibility/falsification workflows.

RGS–PSU–QLG has complementary components: discrete Q6/PRISM geometry and barycentric operators; reconstructed vierbein and teleparallel/torsion diagnostics; spinor, U(1), SU(3), QHD and reduced MHD modules; global and gauged Q-ball boundary-value solvers.

The collaboration should not merge the theories by assertion. Build a shared validation layer: a specified spin/torsion source from RGS can be passed into TorsionLambda's cosmological screening pipeline, while TorsionLambda's consistency gates can be applied to RGS-derived effective stress-energy models.

## 2. Essential distinction: Einstein–Cartan versus propagating torsion

In minimal Einstein–Cartan theory, with the usual Einstein–Cartan–Dirac action and no independent torsion kinetic terms, the Cartan equation is algebraic in torsion: torsion is sourced by the matter spin current. Eliminating algebraic torsion produces an effective local spin–spin interaction (for Dirac matter, an axial-current four-fermion term). It does not automatically introduce an independent propagating torsion wave.

Consequences:
1. Distinguish algebraic Einstein–Cartan torsion from propagating torsion in more general Poincare-gauge or quadratic-curvature/torsion theories.
2. Do not represent contorsion as an unspecified additive tensor in the Einstein equation without deriving the chosen variational formulation and index contractions.
3. A cosmological correction such as H^2=(8 pi G/3)(rho + beta rho_spin) is a phenomenological ansatz until beta and rho_spin are derived from a covariant action, matter spin state, and averaging prescription.
4. If RGS uses a teleparallel formulation, state explicitly whether torsion is the TEGR torsion equivalent to GR or an additional dynamical degree of freedom. These are not interchangeable with Einstein–Cartan torsion.

## 3. Audit items to resolve jointly

### 3.1 Torsion parameter beta

TorsionLambda describes beta as a dimensionless torsion coupling and gives an exploratory interval [-0.5,0.5], associating the bound with causality/hyperbolicity and Holst (1996).

Before using that interval as a physical bound:
- identify the exact action and definition/dimensions of beta;
- derive the principal symbol of the perturbation equations;
- calculate characteristic speeds and hyperbolicity conditions;
- separate exploratory scan priors from derived theoretical bounds and observational posteriors.

Holst's paper concerns the generalized Hilbert–Palatini/Barbero formulation; it should not by itself be cited as establishing a universal numerical bound |beta|<=0.5.

### 3.2 Energy-condition implementation

For a perfect fluid the standard pointwise conditions include:
- WEC: rho >= 0 and rho+p >= 0;
- NEC: rho+p >= 0;
- SEC: rho+p >= 0 and rho+3p >= 0;
- DEC: rho >= 0 and |p| <= rho.

The current WEC implementation checks rho+p but should also check rho>=0. The DEC inequality is not equivalent, in general, to a sound-speed bound c_s^2=dp/drho<=1. NEC violation does not by itself establish closed timelike curves or a generic causality violation. Energy conditions are useful diagnostics, not complete tests of dynamical stability or hyperbolicity.

Recommended joint action: add unit tests for negative rho with positive rho+p, negative pressure saturating DEC, anisotropic stress, and effective spin contributions. Keep local matter conditions distinct from conditions on an effective stress tensor after eliminating torsion.

### 3.3 Spin current and dimensions

The cosmological source rho_spin needs an operational definition. For Dirac matter, start from the canonical spin current and its axial-current reduction in the chosen convention. Record tensor definition and index placement, normalization and units, coarse-graining over an isotropic ensemble, equation of state and conservation law, and relation to the algebraic Cartan torsion solution.

Do not identify a spin density with an energy density without the required coupling scale and dimensions.

### 3.4 Background and perturbation scope

The public documentation describes background FLRW observables and lists full scalar/vector/tensor perturbation equations and CAMB/CLASS likelihood integration as unfinished. Thus the toolkit is presently a pre-simulation screening layer, not a complete Boltzmann solver or observational constraint pipeline. Its background ansatz should be compared with an independently derived EC cosmology before parameter inference is claimed.

## 4. Proposed shared software interface

Create a small adapter contract, not a source-code copy:

    RGSSpinTorsionState:
      scale_factor a
      Hubble_rate H
      energy_density rho
      isotropic_pressure p
      spin_current S[mu,nu,rho]
      axial_current J5[mu]
      torsion T[lambda,mu,nu]
      contorsion K[lambda,mu,nu]
      convention_id
      units_id
      derivation_id

The adapter must reject inputs without convention and units metadata. It should provide a mapping from the selected action to the EC Cartan source; a clearly named phenomenological beta mapping only when independently derived; effective rho/p only after torsion variables have been eliminated or solved; energy-condition diagnostics with the exact tensor and frame stated; limiting-case tests and machine-readable provenance.

No direct coupling of a Q-ball U(1) Noether charge to cosmological spin density is assumed. A bridge requires a specified matter spin current (for example from Dirac fields) and a derivation of its coarse-grained cosmological contribution.

## 5. Joint validation programme

A. **Conventions and action:** freeze metric signature, curvature convention, tetrad orientation, spin-current definition, units, and action. Derive Einstein and Cartan equations variationally.

B. **Algebraic Einstein–Cartan benchmark:** implement the minimal EC-Dirac limit, eliminate torsion analytically, and verify the resulting spin–spin term and GR limit. Compare independent symbolic and numerical evaluations.

C. **Cosmological reduction:** derive the homogeneous/isotropic effective fluid from the specified matter state. Recover standard GR/LambdaCDM when spin coupling is removed. Check conservation/Bianchi identities and dimensional consistency.

D. **Screening and falsification:** run energy-condition tests, parameter scans, limiting cases, and pre-registered fail gates. Correct the WEC implementation and avoid interpreting energy conditions as a substitute for stability.

E. **RGS geometry:** compare the covariant EC torsion solution with RGS vierbein/teleparallel quantities only after identifying the map between the two connections. Do not equate TEGR torsion with EC torsion by notation alone.

F. **Compact objects:** only after A–E, investigate whether Q-ball/superfluid matter with a derived spin current can source torsion in a gravitating compact-object model. A black-to-white-hole transition remains a separate dynamical hypothesis requiring a regular solution and horizon diagnostics.

## 6. Licensing and attribution

The external repository states LGPL-3.0 with a commercial exception. Its checked-in LICENSE appears to contain shortened/simplified license text and a separate commercial-exception section. Before copying, vendoring, or linking source code, obtain clarification from the author about the exact applicable license text, exception, attribution, and contribution terms. The safe initial collaboration is specification-level interoperability and independent reimplementation, with attribution and no copied source.

## 7. Evidence classification

- **Verified from repository contents:** TorsionLambda is a TypeScript pre-simulation screening toolkit with a phenomenological FLRW torsion correction, energy-condition checks, and documented missing perturbation/likelihood stages.
- **Established theoretical background:** minimal Einstein–Cartan torsion is algebraically sourced by spin in the usual action without torsion kinetic terms; integrating it out yields effective spin–spin interactions.
- **Needs correction or derivation:** universal numerical beta bound, spin-density normalization, several energy-condition interpretations, and covariant derivation of the implemented cosmological correction.
- **Open:** a common RGS–PSU–TorsionLambda action, a derived cosmological spin source, propagating torsion sector (if intended), and any compact-object bounce prediction.

## 8. Immediate collaboration deliverable

The first joint milestone should be a short equation-and-convention review, followed by a tiny cross-validation dataset containing (i) GR limit, (ii) algebraic EC spin source, (iii) effective rho/p, (iv) WEC/NEC/SEC/DEC outcomes, and (v) numerical residuals. This is a concrete and falsifiable starting point for collaboration.
