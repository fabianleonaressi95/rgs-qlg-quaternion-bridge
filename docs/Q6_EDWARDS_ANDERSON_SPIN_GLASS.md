# Edwards–Anderson spin glass on Q6

## Scope and evidence status

This module adds a quenched-disorder Ising model to the canonical six-dimensional
hypercube graph Q6. It is deliberately separate from the SU(2) connection and
spinor-current modules.

- **DEFINED:** finite-graph Hamiltonian and quenched Gaussian/bimodal bond ensembles.
- **NUMERICALLY TESTABLE:** gauge invariance, plaquette frustration, finite-size
  replica overlaps, energy and Monte Carlo observables.
- **OPEN:** equilibration and finite-size scaling across a family of graph sizes.
- **NOT CLAIMED:** a thermodynamic spin-glass phase, replica-symmetry breaking,
  equivalence with a three-dimensional experimental spin glass, or a derivation
  from the RGS–PSU vierbein/torsion sector.

## Hamiltonian and quenched disorder

For Ising variables s_i in {-1,+1} on the 64 vertices and independent
couplings J_ij on the 192 edges,

    H_J(s) = - sum_(ij) J_ij s_i s_j.

Supported ensembles are Gaussian J_ij ~ Normal(0,J0^2) and bimodal
J_ij = +/- J0 with equal probability. A seed selects one disorder realization.
The couplings remain fixed throughout thermal sampling: they are quenched, not
redrawn at each sweep.

## Frustration and gauge invariance

For each of the 240 oriented square faces,

    F_p = product_(ij in boundary(p)) sign(J_ij) in {-1,+1}.

F_p=-1 marks a frustrated square. Under the local Z2 transformation

    s_i -> eps_i s_i,       J_ij -> eps_i J_ij eps_j,

the Hamiltonian and every F_p are invariant. This is a Z2 bond-gauge invariance
of the Ising model; it is not the SU(2) gauge symmetry of the existing Wilson
link construction.

## Replica observables

Two independent parallel-tempering ladders use the same frozen disorder sample.
At each temperature, synchronized measurements define

    q_12 = (1/64) sum_i s_i^(1) s_i^(2),
    chi_SG = 64 <q_12^2>,

and the Binder cumulant

    U_4 = 1 - <q_12^4> / (3 <q_12^2>^2).

These are finite-size estimators. Reliable conclusions require multiple
independent disorder samples, autocorrelation/equilibration checks, temperature
round trips, and scaling over a sequence of graph sizes. A single Q6 sample or
a single Binder curve cannot establish a phase transition.

## Exact and approximate optimization

The exact_ground_state function uses Gray-code enumeration and is intended only
for small induced graphs. It refuses the full 64-vertex graph by default because
2^64 configurations are not a practical exact audit. The full Q6 dynamics uses
parallel tempering; its lowest sampled energy is not certified as the ground
state energy.

## Reproduction

    PYTHONPATH=. python -m pytest tests/test_rgs_q6_spin_glass.py -q

The API is Q6SpinGlass.random, model.frustration, model.gauge_transform,
exact_ground_state and parallel_tempering in src/rgs_q6_spin_glass.py.

## Literature starting points

- S. F. Edwards and P. W. Anderson, “Theory of spin glasses,” J. Phys. F 5
  (1975), 965–974, DOI: 10.1088/0305-4608/5/5/017.
- D. Sherrington and S. Kirkpatrick, “Solvable Model of a Spin-Glass,”
  Phys. Rev. Lett. 35 (1975), 1792–1796, DOI: 10.1103/PhysRevLett.35.1792.
- G. Parisi, “Infinite Number of Order Parameters for Spin-Glasses,”
  Phys. Rev. Lett. 43 (1979), 1754–1757, DOI: 10.1103/PhysRevLett.43.1754.
- K. Binder and A. P. Young, “Spin glasses: Experimental facts, theoretical
  concepts, and open questions,” Rev. Mod. Phys. 58 (1986), 801–976,
  DOI: 10.1103/PhysRevLett.58.801.
- M. Mézard, “Spin glass theory and its new challenge: structured disorder,”
  Indian J. Phys. (2024), DOI: 10.1007/s12648-023-03029-8.

The references motivate the model and diagnostics; they do not validate a
physical identification of Q6 with a material spin glass.
