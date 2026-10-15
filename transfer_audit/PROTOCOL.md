# Structure-matched small-molecule transfer audit

Protocol ID: `structure_matched_transfer_exploratory_2026_10_05`.
Protocol amendment: 2026-10-05 (Asia/Shanghai), issued before computing any new-cohort rankings.

Status: **post hoc exploratory transfer audit**, not preregistered or confirmatory. The original 14-formula benchmark and 49 small-molecule tasks have already been analysed. A first attempted 60-task selection on 2026-10-04 was found to omit required approval-year and ATC-stratum eligibility rules. Ranking was not run. Its invalidated runtime files are retained under `runtime_invalidated_60task_selection/` for provenance only: five tasks have aggregate activity and PubMed response files, and one additional task has an activity response but no PubMed aggregate; target-metadata completeness was not certified. These data are excluded from all analyses. The corrected cohort below reuses the original eligibility criteria and is selected without reference to ranking outputs or external target labels.

## Question

Can the unchanged evidence-gated workflow be applied to all remaining eligible small-molecule tasks with exact structure identity resolution, and what are the rankability, missingness, and ranking-change distributions? A secondary, descriptive DrugCentral FDA-label audit may assess whether externally sourced historical annotations can be mapped to the ranked universe; it is not a gold standard and may share upstream evidence with ChEMBL or PubMed.

## Frozen frame and cohort selection

- Source frame: the exact 4,225 ChEMBL molecule records already frozen for the original transfer audit; preserve its manifest and page hashes.
- Reapply the original eligibility without change: `max_phase=4`; `molecule_type=Small molecule`; `structure_type=MOL`; first approval 2010–2022 inclusive; self-parent molecule; nonempty canonical SMILES; and ATC stratum assigned by priority C, then J, then N.
- Identity: require exact Standard InChIKey equality between a ChEMBL parent and DrugCentral `structures.smiles.tsv` (2021-09-01 source path); the key must map to exactly one ChEMBL parent and one DrugCentral structure ID. Do not use names to resolve ambiguity.
- Exclude all 49 prior ChEMBL task IDs. Include **every** remaining eligible task with an unambiguous structure link; do not cap the cohort, sample according to activity/label availability, or replace tasks. Retain source-frame order within C/J/N strata using the original hash seed `EthnoEvidence-transfer-2026-10-01|` and report all strata and shortfalls.
- External target labels are examined only after this ID list is frozen. For the secondary audit, retain unique DrugCentral rows with `ORGANISM=Homo sapiens`, `MOA=1`, `MOA_SOURCE=DRUG LABEL`, an FDA label URL year from 2003–2022, and matching DrugCentral structure ID. Post-cutoff labels are ineligible; no label means missing reference data, never a negative.

## Workflow and endpoints

Use the existing small-molecule adapter and ranking specification unchanged: cutoff 2022-12-31, primary pChEMBL threshold 5.0, original measured/literature eligibility, fixed candidate universe, and the evidence-gated, measured-only, and ungated-union comparators. Preserve retrieval URLs, timestamps, page counts, source hashes, errors, and incomplete states. Do not tune weights or thresholds against this cohort.

Primary operational summaries: selected task count; retrieval-complete tasks by source; rankable tasks; candidate-universe size; evidence/source counts; Top-10 replacement and Jaccard distributions across the fixed-universe rankers; threshold sensitivity; and all missingness/failure reasons. These are workflow-transfer and decision-change descriptors, not accuracy metrics.

Secondary external-label summaries, only when identity/date checks pass: label coverage, number of dated target pairs in the candidate universe and Top-10, and record-level source overlap. Do not report sensitivity, recall, precision, AUC, NDCG, superiority, or biological validation without a separately designed and independently adjudicated reference set.

## Limitations and release boundary

The cohort is post hoc relative to the original project. The DrugCentral structure map is dated 2021-09-01; ChEMBL and PubMed are current database contents filtered retrospectively by publication/document year rather than reconstructed historical snapshots. FDA label annotations may overlap upstream evidence. DrugCentral CC BY-SA raw records stay outside the MIT public release unless license compatibility is separately resolved. Preserve all null results, exclusions, and deviations; do not pool this exploratory cohort with the frozen formula benchmark as if they were exchangeable observations.

## Source provenance

The DrugCentral download page describes the SMILES/InChI file as including structure identifiers and DrugCentral IDs; the local 2021-09-01 structures file was retrieved from `https://unmtid-dbs.net/download/DrugCentral/2021_09_01/structures.smiles.tsv`. The interactions TSV is described as incorporating literature, drug-label, and external data-source target interactions. DrugCentral's official license page links to CC BY-SA 4.0. The two local input SHA-256 values are recorded in `cohort.json` and `cohort_and_reference_audit.json`; licensed source rows are retained only in this local audit workspace, not copied into the MIT release candidate.
