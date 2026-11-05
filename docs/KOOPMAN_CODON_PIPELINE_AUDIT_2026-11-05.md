# Koopman codon-pipeline audit (2026-11-05)

## Findings

1. experiments/quantum_models/koopman2.py previously loaded base.py with exec(open(...).read()). This executes source text and resolves the file relative to the caller working directory. It now imports the sibling module normally.
2. base.py previously opened codon_train.csv relative to the caller working directory. The loader now resolves the file in this order: explicit CODON_TRAIN_CSV, next to base.py, then current directory. An explicit but invalid environment path fails rather than silently selecting another file.
3. The repository tree inspected for this audit contains no codon_train.csv. Consequently, this PR cannot execute the biological analysis or validate its schema/results. The dataset must be supplied separately; no replacement dataset or digest is fabricated.
4. The comment in base.py says “split by organism”, but the code performs rng.permutation(len(X)) followed by an 80/20 row split. This is a random row split, not an organism-held-out split. The visible code exposes Kingdom, not a confirmed organism identifier. Do not interpret this as out-of-organism generalization. A valid grouped split requires a documented organism/taxon ID column and a group-level split implementation.
5. The null comparisons in koopman2.py reuse the fixed held-out pair set Pte for evaluation while changing training pairs. This is a fixed-test design only if stated explicitly; the nulls do not independently generate test sets. This audit does not upgrade these results to inferential evidence.

## Evidence classification

- Verified from source: dynamic exec, CWD-dependent CSV loading, random row split, fixed held-out set reused in null evaluations.
- Not verified: dataset identity, row-to-organism mapping, data integrity, scientific conclusions, numerical outputs, and organism-level generalization.
- Required next input: exact codon_train.csv (or a stable source URL plus checksum) and its data dictionary, including a unique organism/taxon ID.