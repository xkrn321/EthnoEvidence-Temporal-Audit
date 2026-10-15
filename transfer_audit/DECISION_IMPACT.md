# Task-level decision-impact audit

Protocol: `structure_matched_transfer_exploratory_2026_10_05`  
Source: `task_level_audit.tsv` and `runtime/rankings.json`  
Status: post hoc exploratory description; not an accuracy evaluation.

## Top-10 changes under unchanged ranking rules

| Contrast | Rankable tasks with ≥1 replacement | Replacement range among changed tasks | Most changed task |
|---|---:|---:|---|
| Evidence-gated vs ungated union | 6/36 (16.7%) | 1–4 | CHEMBL231813: 4/10 replaced; Jaccard 0.429 |
| Evidence-gated vs measured-only | 9/36 (25.0%) | 1–9 | CHEMBL40: 9/10 replaced; Jaccard 0.053 |

The median replacement count across all 36 tasks was zero for both contrasts. This is a heterogeneous response: generally stable Top-10 membership, with larger changes in a small number of tasks. Do not describe it as a general or frequent ranking improvement.

## Stratification by measured-data availability

| Rankable task stratum | n | Gated vs ungated: tasks changed | Gated vs measured-only: tasks changed | Interpretation |
|---|---:|---:|---:|---|
| At least one eligible measured-supported target | 24 | 6/24 (25.0%) | 7/24 (29.2%) | Both fixed-universe contrasts are calculable; median replacements remain zero. |
| No eligible measured-supported target | 12 | 0/12 | 2/12 (16.7%) | Measured-only contains no measured scores and falls back to alphabetical ordering; treat its changes as a tie-break artifact, not an informative comparator result. |

These are descriptive post hoc strata, not confirmatory subgroups. The more interpretable decision-impact result is the gated-versus-ungated change among tasks that do have measured support; it still says only that the rule changes a minority of lists, not that it improves correctness.

## Interpretable extreme case

CHEMBL40 had 79 candidates, zero eligible measured-supported targets, and 79 literature-supported targets. The “measured-only” comparator was computed on the shared fixed candidate universe; because all measured scores were absent, its order fell back to alphabetical tie-breaking. Its Top 10 shared one target with the gated Top 10 (9 replacements; Jaccard 0.053). The gated-only entries were CYP1A1, CYP1A2, CYP2B6, CYP2E1, PCNA, CYP2C19, CYP3A4, CYP2C9, and UGT1A; the comparator-only entries were ABCB11, ABCC2, ABCC4, ABCC6, ABCG2, ALDH1A3, ALDH3A1, ARAF, and BDNF. This is an auditable response under the declared rules, but the comparator is uninformative for this task and the contrast does not identify which Top 10 is biologically correct. No FDA-label target pair was available for this task.

Another case, CHEMBL290916, had 28 candidates, no eligible measured-supported targets, and 28 literature-supported targets; four Top-10 entries changed relative to measured-only (Jaccard 0.429). These examples should be presented as illustrative, not as a selected estimate of average performance.

## Candidate generation versus ranking

The separate DrugCentral label audit found seven of 21 historical label pairs outside the candidate universe. Reranking cannot recover a target that candidate generation never supplies. A future benchmark should therefore measure candidate-universe coverage separately from within-universe ranking, using a reference set that is both independent of inputs and adjudicated without access to rankings.

## Reproducibility

The transfer summary script rebuilt 76,958 evidence objects from stored raw responses and exactly reconstructed all feature objects, ranking lists, and task-level decision audits. Removing historical FDA-label targets from task metadata left every ranking, feature, and decision audit unchanged, verifying that reference labels are not passed into the implemented ranking calculation. This does not rule out indirect source dependence or selection bias outside that calculation. `top10_membership_changes.tsv` lists every added and removed member for all 15 task–comparator contrasts that changed. This is a verified computational reproducibility result for the saved cohort and inputs, not independent replication or biological validation. Core transfer, ranking, and DrugCentral-reference tests passed (28 tests). A first repository-wide pytest invocation from the wrong import context failed collection; it is not counted as a test failure of the code and the full suite was not re-run in its configured root context during this audit.
