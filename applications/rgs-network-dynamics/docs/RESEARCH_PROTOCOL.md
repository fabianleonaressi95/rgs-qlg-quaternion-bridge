# Competence-Based Self-Management — research protocol

**Version:** 1.0  
**Date:** 2026-11-19  
**Status:** research protocol; empirical validation not yet performed.

## 1. Research question

Are organizations structured around distributed competencies more resilient, adaptable or productive than organizations structured primarily around roles and functions?

The study does not assume that self-management is superior to hierarchy. It compares organizational mechanisms under explicit, measurable conditions.

## 2. Cases and analytical roles

| Case | Unit of analysis | Mechanism | Evidence status at project start |
|---|---|---|---|
| Practices / Sun Microsystems Italia | Competence Practice | Shared expertise and focal point | Primarily testimonial; independent historical data missing |
| Buurtzorg | Self-managing care team | Professional autonomy, coaches and support infrastructure | Published case research and institutional reports |
| Morning Star | Colleague / commitment | Peer agreements and individual responsibility | Harvard teaching case and secondary analysis |
| FAVI | Customer-oriented production unit | Autonomous mini-factories | Field research and institutional case |
| Sun Hydraulics | Operating unit | Reduced formal hierarchy | Historical teaching case |
| Zappos / Holacracy | Circle and role | Formalized distributed governance | Published case analyses; transition difficulties documented |
| Teal organizations | Organization-level interpretive category | Self-management, wholeness, evolutionary purpose | Theoretical framework; cases are heterogeneous |

These cases are not interchangeable. A Practice, a nursing team, a production cell and a Holacracy circle are different units of analysis. Comparisons must record organizational level and period.

## 3. Hypotheses

- **H1 — Competence redundancy:** greater distribution of critical competencies is associated with lower capability loss after an agent becomes unavailable.
- **H2 — Competence concentration:** higher concentration of a critical skill is associated with greater operational vulnerability.
- **H3 — Knowledge transfer:** stronger documented knowledge-transfer relations are associated with lower knowledge loss over time.
- **H4 — Decision autonomy:** greater autonomy may reduce decision latency when tasks are modular and responsibilities are clear.
- **H5 — Coordination:** removing hierarchical levels does not necessarily reduce coordination costs; the effect depends on task complexity and decision rules.
- **H6 — Structural resilience:** resilience depends on competence distribution and alternative collaboration paths, not simply on headcount or graph connectivity.
- **H7 — Autonomy × competence:** the effect of autonomy on outcomes depends on whether agents have the skills and information needed to decide.

All are hypotheses, not established findings. Sector, size, time period, task complexity and workforce composition are potential confounders.

## 4. Multilayer representation

Represent an organizational unit by \(\mathcal G=(G_C,G_W,G_K,G_A)\):

- \(G_C\): agent-to-competence relation;
- \(G_W\): operational collaboration;
- \(G_K\): directed knowledge transfer;
- \(G_A\): formal authority, approval and delegation.

Do not merge these layers without an explicit model. A focal point is not coded as a hierarchical superior unless the evidence documents formal authority.

## 5. Competence matrix

Let \(C\in[0,1]^{n\times m}\), where \(C_{ik}\) is the documented proficiency of agent \(i\) in skill \(k\). Use a rubric with observable anchors; do not infer proficiency from job title alone.

For threshold \(\tau\), the number of qualified agents for skill \(k\) is
\[
N_k(\tau)=\sum_i\mathbf 1(C_{ik}\geq\tau).
\]
For \(s_k=\sum_i C_{ik}>0\), define
\[
D_{ik}=C_{ik}/s_k,\qquad H_k=\sum_iD_{ik}^2.
\]
\(H_k\) is a Herfindahl-style concentration index. It is undefined when no agent has the skill. Its reciprocal is an effective-number interpretation, not a direct measure of substitutability.

A skill is not operationally redundant merely because several people list it: verification must establish that they can independently perform the relevant task.

## 6. Knowledge dynamics (exploratory model)

Let \(X_t\in[0,1]^{n\times m}\) be the agent-by-skill knowledge state and \(P_K\) a row-stochastic knowledge-mixing matrix. An exploratory discrete model is
\[
X_{t+1}=\Pi_{[0,1]}\left[X_t+
\eta(P_KX_t-X_t)-\delta X_t+U_t\right].
\]
Here \(\eta\) is a mixing coefficient, \(\delta\) a decay parameter and \(U_t\) an external learning input. This is a phenomenological modelling choice, not an established law of organizational learning. Parameters must be estimated from longitudinal observations or explored transparently through sensitivity analysis.

## 7. Outcomes and indicators

Record separately:
- decision autonomy: share of decisions not requiring superior approval;
- competence coverage: share of critical skills covered above a declared proficiency threshold;
- competence concentration: \(H_k\) by skill;
- knowledge transfer: observed directed transfer events, not merely declared ties;
- coordination: elapsed time and number of handoffs per task;
- resilience: task capability retained after simulated or observed unavailability;
- performance: quality, cost, productivity and satisfaction, with operational definitions.

Graph connectivity is not business continuity. Knowledge-state diffusion is not output quality. No composite score should be introduced unless weights are preregistered and justified.

## 8. Resilience experiment

Define a task by the set of skills required to execute it. A task remains feasible if every required skill is available among the remaining agents above threshold \(\tau\). Compare:
1. random agent removal;
2. removal of high-degree agents;
3. removal of agents with unique or highly concentrated critical skills.

Report baseline coverage and post-removal coverage separately. If baseline coverage is zero, the retention ratio is undefined and must not be silently set to zero.

## 9. Mixed-method design

For each case collect:
- period, sector, organization size and unit of analysis;
- org charts and role descriptions;
- decision and delegation rules;
- competence matrices and training records;
- collaboration and knowledge-transfer records;
- task resolution times and handoffs;
- quality, cost, productivity and workforce outcomes;
- interviews with workers, coordinators and managers.

Interview prompts:
1. Who decides when a problem occurs?
2. Who has the necessary competence?
3. How is the right person found?
4. What happens when that person is unavailable?
5. How is the knowledge retained after the problem is solved?

Use common coding rules and preserve contradictory accounts.

## 10. Evidence grading

- **A — primary verifiable:** contemporaneous company records, administrative data, auditable operational records.
- **B — independent research:** peer-reviewed studies or institutional research with described methods.
- **C — testimonial:** interviews and participant accounts.
- **D — interpretive:** theoretical reconstructions or descriptions without independent verification.

Economic or causal claims cannot be certified from C/D sources alone. Synthetic simulations are never empirical validation.

## 11. Practices-specific evidence gate

Prioritize original 1998-era documents and independent corroboration:
- date and context of introduction at Sun Microsystems Italia;
- number and scope of Practices;
- membership and focal-point selection rules;
- actual decision rights of focal points;
- routing and assignment of requests;
- substitution and knowledge-retention procedures;
- training and international rollout records;
- independent accounts from former employees.

The claims of adoption in 30 versus more than 40 countries and a 20% profit advantage remain unverified until source documents, periods, comparison groups and accounting definitions are supplied. No causal effect is inferred from the testimonial claim.

## 12. Statistical analysis and limitations

Start with descriptive statistics, missingness analysis, uncertainty intervals and sensitivity checks. Longitudinal regression or causal designs are appropriate only when the number and structure of observations support them. Adjust for sector, scale, period, task complexity and workforce composition.

The initial case set is small and heterogeneous. It supports structured comparison and hypothesis generation, not universal generalization.

## 13. Completion gates

1. **Documentary:** every central case claim has a traceable source and evidence grade.
2. **Network:** node/edge definitions and data provenance are explicit.
3. **Mathematical:** edge cases, normalization, thresholds and removal algorithms are tested.
4. **Empirical:** indicators are checked against observed operations.
5. **Comparative:** hypotheses are evaluated with confounders and uncertainty reported.

No claim of Practices superiority or RGS-specific modelling advantage is permitted before the relevant empirical or computational gates pass.
