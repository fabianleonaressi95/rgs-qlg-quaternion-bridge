# RGS–PSU core parse repair and finite 600-cell moment audit

Date: 2026-11-04

## Changes

- Removed a duplicated, stale article-generation block from `RGS_PSU_closure.py` that began after the first complete report had already been written. The duplicate block caused an unexpected-indentation syntax error and contained conflicting claims about a converged BVP and the teleparallel identity.
- Replaced the malformed equidistribution script with an explicit 120-vertex 600-cell coordinate construction and a reproducible finite-moment audit.
- Added tests for parsing, vertex count, unit norms, moment residuals and the quotient distinction.

## Mathematical correction

A finite orbit of the binary icosahedral group is not the continuous Haar measure on (S^3). Symmetry or agreement of a finite number of moments does not establish weak equidistribution of a sequence of refined point measures.

The 600-cell vertex set contains 120 points. It is the orbit of the identity under the binary icosahedral group action. Under the left quotient (S^3/2I), the points of that single orbit all represent one quotient point. This must not be confused with the antipodal quotient (S^3/{+1,-1}), which is (mathbb{RP}^3).

The new script verifies numerically that the 120 coordinates:
- lie on the unit (S^3);
- have vanishing first moment;
- match the normalized (S^3) second and fourth Haar moments to floating-point tolerance.

These are finite moment checks (consistent with a spherical-design property). They do **not** establish continuous equidistribution, refinement convergence, Hausdorff convergence, or a Mosco limit.

## Evidence classification

- **NUMERICALLY VERIFIED:** coordinate count, uniqueness, unit norms and first/second/fourth moment residuals for this finite 120-point configuration.
- **PROVED by construction:** the generated coordinate list has 120 unique entries, assuming the exact algebraic coordinate prescription implemented.
- **OPEN:** weak convergence of refined measures, density/Hausdorff convergence, quotient-space refinement dynamics and continuum operator convergence.
- **Not claimed:** a proof of Weyl equidistribution for the finite 2I orbit.

The closure script repair is a syntax/duplicate-block correction only. It does not complete the teleparallel action, resolve tetrad/spin-connection conventions, derive the NLED sector or validate a coupled wormhole BVP.
