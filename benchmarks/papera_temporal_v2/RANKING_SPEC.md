# Temporal resource v2 ranking specification

## Candidate universe

For every formula, the primary candidate universe is the union of approved
human gene symbols supported by either (a) document-year-constrained ChEMBL
activity records dated no later than 2022 or (b) high-specificity pre-cutoff
PubMed mentions. Current SwissTargetPrediction (STP) results are retained only
for a non-temporal comparator; they cannot create a primary candidate.

## Evidence-gated ordering

Primary candidates are ordered lexicographically by:

1. strongest eligible temporal tier: strong measured, weak measured, then
   pre-cutoff literature-only;
2. number of distinct recorded source classes among measured evidence and
   pre-cutoff literature;
3. number of distinct formula marker compounds supporting measured evidence;
4. maximum retained pChEMBL value; and
5. number of pre-cutoff article mentions, then alphabetical HGNC symbol.

No current STP signal changes this primary key.

## Comparators and ablations

- `pre_cutoff_literature_frequency`: pre-cutoff mention frequency only.
- `cutoff_chembl_measured_only`: pure measured-evidence tier, pChEMBL, and
  measured compound coverage; STP and literature cannot break ties.
- `stp_only_non_temporal_comparator`: current STP probability on the same
  primary candidate universe; it does not create candidates.
- `ungated_union`: equal-weight temporal measured strength, literature
  frequency, and measured compound coverage, with missing components set to 0.
- `seeded_random`: deterministic shuffle of the same universe.
- `without_measured_tier`, `without_pre_cutoff_literature`, and
  `without_missingness_penalty`: distinct ablations only.

No duplicate alias of `ungated_union` is reported as an independent ablation.
