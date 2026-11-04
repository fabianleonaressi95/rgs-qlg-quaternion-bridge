# Gödel–RGS–Wyler: an ontological-to-geometric research bridge

Status: research programme / conjectural bridge. This note does **not** claim a
derivation of the fine-structure constant or a physical proof of Gödel's
ontological conclusion.

## 1. Three different mathematical objects

### 1.1 Gödel: modal admissibility of properties

In a representative formulation, Gödel's proof uses a predicate \(P(\varphi)\)
for positive properties, quantified modal logic, an essence relation, and
necessary existence. A simplified axiom pattern includes:

- positivity is exclusive with positivity of the negation;
- positivity is closed under necessary entailment;
- positivity is necessarily positive;
- necessary existence is positive.

The conclusion follows relative to the chosen axioms and modal semantics.
It is not a measure on a space of physical configurations. In particular,
"positive" in Gödel's technical predicate must not be silently identified
with positive energy, positive-definite action, probability, or an RGS
admissibility score.

### 1.2 RGS: a space of admissible structures

Let \(\mathfrak X\) be a specified configuration space (for example, finite
Q6 cochain/gauge data at a fixed resolution), with a group \(G\) acting on it.
Let \(\mathcal C=\{C_a\}\) be explicit structural constraints, such as
incidence identities, equivariance, unitarity, and chosen boundary data.
Define
\[
  \mathfrak X_{\rm adm}=\{x\in\mathfrak X:C_a(x)=0\ \forall a\}/G.
\]
This is a concrete quotient only after the configuration space, constraints,
and group action have been specified. It is not automatically the set of
Gödel-positive properties or possible worlds.

A property \(p\) can be represented by a predicate \(p:\mathfrak X\to\{0,1\}\).
A useful RGS analogue of "joint consistency" is the non-emptiness test
\[
  \bigcap_{p\in\mathcal P_+}\{x:p(x)=1\}\cap\mathfrak X_{\rm adm}
  \ne\varnothing.
\]
This is a proposed bridge condition, not one of Gödel's theorems. It is
stronger and more operational than merely assigning a label "positive".

### 1.3 Wyler: invariant geometry and normalization

Wyler's proposal expresses a numerical approximation to \(\alpha\) through
geometric quantities associated with bounded symmetric domains. The
mathematical challenge is not just evaluating volumes: one must justify the
domain, invariant measure, normalization, scale/radius, and the physical
identification with a gauge coupling.

Robertson's 1971 critique specifically observed that the agreement depended
on setting a radius to one without a known reason; Gilmore's 1972 response
argued for the distinguished normalization in terms of coset representatives.
This exchange means the normalization issue must be treated explicitly,
not hidden inside a numerical match.

## 2. A disciplined bridge

The proposed chain is

\[
\boxed{
\text{modal axioms}
\;\longrightarrow\;
\text{consistency/admissibility}
\;\longrightarrow\;
\text{quotient geometry}
\;\longrightarrow\;
\text{invariant measure}
\;\longrightarrow\;
\text{effective gauge action}
}
\]

Each arrow requires a separate theorem or explicit modelling assumption.

### Step A — modal layer

Choose a Kripke frame \((W,R)\), a domain policy, and an interpretation of
properties. State whether the intended logic is S5 or another normal modal
system. Check consistency of the positivity axioms in the chosen semantics.
Do not infer physical existence from formal necessitation.

### Step B — structural layer

For each resolution \(n\), define \(\mathfrak X_n\), its gauge/symmetry group
\(G_n\), and constraints \(C_{n,a}\). In the Q6 programme these may include
\(d^2=0\), icosahedral equivariance, SU(2) link covariance, and Clifford
relations, but these constraints describe mathematical structures, not
ontological positivity.

A possible semantics map is
\[
  \llbracket p\rrbracket_n
  =\{[x]\in\mathfrak X_{n,\rm adm}:p(x)=1\}.
\]
A family of properties is jointly realizable precisely when the intersection
of their extensions is nonempty. This gives an operational analogue of
compatibility, not Gödel's original predicate \(P\).

### Step C — geometric layer

If the quotient admits a suitable manifold/orbifold or a well-defined
stratified measure space, specify its metric and invariant volume form.
For a homogeneous space \(G/K\), identify the group, isotropy subgroup,
domain, and normalization convention. Prove that the measure is invariant
and finite (or specify a renormalized measure if it is not).

A candidate partition measure is
\[
 d\mu_n([x])=Z_n^{-1}e^{-S_n([x])}\,d\nu_n([x]),
\]
where \(d\nu_n\) is independently defined. The action \(S_n\) must not be
chosen after seeing the desired value of \(\alpha\).

### Step D — physical matching

Integrate out the specified degrees of freedom and demonstrate a controlled
continuum/low-energy limit:
\[
 \Gamma[A]\supset
 \frac{1}{4g^2(\mu)}
 \int d^4x\,\sqrt{|g|}\,F^a_{\mu\nu}F^{a\mu\nu}.
\]
Then specify gauge group, generator normalization, renormalization scheme,
scale \(\mu\), matter content, and matching to the Thomson-limit fine
structure constant. Only after these choices can a geometric coefficient be
compared with \(\alpha\).

## 3. What Q6 currently supplies—and what it does not

The repository now contains finite cubical cochains, an equivariant lift of
the icosahedral action, SU(2) parallel transport on edges, and an exterior
Clifford module. These supply finite-dimensional algebraic ingredients.

They do **not** yet supply:
- a modal semantics or a Gödel-positive-property interpretation;
- a canonical identification of the Q6 quotient with a Wyler bounded
  symmetric domain;
- a canonical invariant measure and its scale normalization;
- a continuum limit yielding a four-dimensional gauge action;
- a renormalized prediction for \(\alpha\).

In particular, the 64-dimensional exterior Clifford module is not by itself
a physical spinor field, and internal SU(2) link transport is not a Lorentz
spin connection.

## 4. Falsifiable proof obligations

1. **Semantics:** give a precise model of the proposed admissibility predicate
   and prove its consistency independently of any target coupling.
2. **Functoriality:** show that admissible maps preserve the relevant
   structural predicates and descend to the quotient.
3. **Geometry:** construct the candidate homogeneous space and prove the
   claimed invariant measure, including all constants and scale choices.
4. **Uniqueness:** show that the domain and normalization are forced by
   structural axioms, rather than selected to fit \(\alpha\).
5. **Dynamics:** derive the coefficient of \(F^2\) from one specified action.
6. **Continuum/renormalization:** prove or numerically control convergence and
   match the physical coupling at a declared scale.
7. **Robustness:** vary admissible discretizations and prove the prediction is
   stable under refinement.

A failed obligation is informative: it rules out that particular bridge.

## 5. Epistemic ledger

| Statement | Status |
|---|---|
| Gödel's modal argument is derivable relative to its axioms and semantics | Established formal result, conditional on assumptions |
| RGS finite Q6 structures can encode explicit admissibility constraints | Implemented finite mathematics |
| Gödel positivity equals RGS admissibility | Not established; presently a category error if asserted directly |
| Q6 quotient is a Wyler domain | Open conjecture; no identification supplied |
| A canonical measure fixes Wyler's normalization | Open |
| The measure produces the gauge kinetic coefficient | Open |
| The resulting coefficient predicts measured \(\alpha\) | Not demonstrated |

## 6. Terminological guardrail: Wyler versus Weyl

This note concerns **Armand Wyler** and the geometric proposal for the
fine-structure constant. **Hermann Weyl** is a separate route involving
symmetry, gauge ideas, and geometry. Weyl's gauge-geometric programme may
provide useful mathematical tools, but it must not be conflated with Wyler's
volume formula.

## References

- K. Gödel, ontological proof notes; see the scholarly reconstruction in
  Jordan Howard Sobel, *Logic and Theism*, Cambridge University Press,
  chapter "Kurt Gödel's Ontologischer Beweis".
- Stanford Encyclopedia of Philosophy, "Ontological Arguments",
  section on Gödel: https://plato.stanford.edu/entries/ontological-arguments/
- B. Robertson, "Wyler's Expression for the Fine-Structure Constant
  \(\alpha\)", *Physical Review Letters* 27 (1971), 1545,
  https://doi.org/10.1103/PhysRevLett.27.1545
- R. Gilmore, "Scaling of Wyler's Expression for \(\alpha\)",
  *Physical Review Letters* 28 (1972), 462,
  https://doi.org/10.1103/PhysRevLett.28.462


## Primary literature and distinction of Gödel results

For Gödel's modal ontological argument, see the historical reconstruction by Kanckos and Lethen (2019), *The Development of Gödel's Ontological Proof*, Review of Symbolic Logic: https://www.cambridge.org/core/journals/review-of-symbolic-logic/article/development-of-godels-ontological-proof/728D16FEF40BA2054055E378C879F369, and Benzmüller and Scott (2025), *Notes on Gödel's and Scott's variants of the ontological argument*, https://doi.org/10.1007/s00605-025-02078-x. For the modal frame semantics, see Blackburn, de Rijke and Venema, *Modal Logic* (Cambridge, 2001).

Gödel's 1931 incompleteness paper is a separate result: *Über formal unentscheidbare Sätze der Principia Mathematica und verwandter Systeme I*, https://doi.org/10.1007/BF01700692. It must not be conflated with the later ontological argument.

For the geometric expression and its radius issue, see Wyler (1969), C. R. Acad. Sci. Paris 269A, 743–745; Robertson (1971), https://doi.org/10.1103/PhysRevLett.27.1545; and Gilmore (1972), https://doi.org/10.1103/PhysRevLett.28.462. These are historical/theoretical sources for the expression and normalization debate, not a validation of the proposed RGS-to-alpha derivation.

The present note remains a research bridge: no cited source proves the arrows from local admissibility to a D5 invariant measure, or from that measure to a renormalized physical gauge coupling.