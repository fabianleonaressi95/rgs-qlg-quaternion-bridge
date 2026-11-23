# Data dictionary

All real organizational data must be authorized, minimized and pseudonymized. Synthetic data must be marked `synthetic` and stored separately from empirical data.

## agents.csv

| Field | Type | Meaning |
|---|---|---|
| agent_id | string | Pseudonymous agent identifier |
| unit_id | string | Team / organizational unit |
| period_id | string | Observation period |
| role_type | string | Role description or normalized role class |
| autonomy_level | float [0,1] | Evidence-based autonomy score |
| active | boolean | Agent available in period |

## skills.csv

| Field | Type | Meaning |
|---|---|---|
| skill_id | string | Stable skill identifier |
| skill_name | string | Skill label |
| criticality | float [0,1] | Declared task-criticality score |
| threshold | float [0,1] | Proficiency threshold for operational coverage |

## agent_skills.csv

| Field | Type | Meaning |
|---|---|---|
| agent_id | string | Agent key |
| skill_id | string | Skill key |
| proficiency | float [0,1] | Assessed proficiency |
| evidence_type | string | Assessment, certification, task record, interview, etc. |
| verified | boolean | Whether the proficiency was independently checked |

## relations.csv

| Field | Type | Meaning |
|---|---|---|
| source_id | string | Source agent |
| target_id | string | Target agent |
| relation_type | string | collaboration, knowledge_transfer, authority, approval |
| weight | float >=0 | Frequency or intensity with declared unit |
| period_id | string | Observation period |
| evidence | string | Provenance reference |

## decisions.csv

| Field | Type | Meaning |
|---|---|---|
| decision_id | string | Decision key |
| unit_id | string | Organizational unit |
| decision_type | string | Decision category |
| decision_owner | string | Pseudonymous owner |
| approval_required | boolean | Whether another approval was required |
| resolution_time | float >=0 | Time to resolution; unit must be declared |

## outcomes.csv

| Field | Type | Meaning |
|---|---|---|
| unit_id | string | Organizational unit |
| period_id | string | Observation period |
| productivity | float or null | Defined output/input measure |
| quality | float or null | Defined quality measure |
| cost | float or null | Cost with currency and period |
| satisfaction | float or null | Survey measure and instrument |
| turnover | float or null | Rate and denominator |
| absenteeism | float or null | Rate and denominator |

## Provenance requirements

Each imported dataset must have a companion manifest recording source, access date, owner, permission basis, collection method, time coverage, transformations, missing-data handling and whether the data are synthetic or empirical. Never store direct identifiers unless strictly necessary and authorized.
