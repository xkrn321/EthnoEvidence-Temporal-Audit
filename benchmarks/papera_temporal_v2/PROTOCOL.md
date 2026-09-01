# Temporal evaluation resource v2

Amended on 29 August 2026 before submission, after an internal data-integrity
audit identified temporal-date inconsistencies and an overbroad automated proxy
in v1. Version 1 remains preserved as a superseded audit artifact and must not
be used for manuscript figures or claims.

## Scope and claim boundary

This is a software and resource study. It documents candidate generation,
evidence provenance, screening aids, ranking reproducibility, and missing-data
handling. It does **not** estimate predictive performance until two blinded
reviewers have completed full-text screening and target adjudication.

## Cohort and temporal boundary

- Cohort: 14 predeclared classical formulas in the versioned registry.
- Excluded prospective deployment case: retained outside this benchmark and
  reserved for separately authorized staged validation.
- Development evidence: canonical first publication year no later than 2022.
- Primary held-out screening window: 1 January 2023 to 31 December 2025.
- Retrieval: one broad PubMed query from 1900 through 2025 per formula, followed
  by one shared canonical-date policy. This prevents ESearch's publication-date
  field from being mixed with a different parsed year.
- Canonical date: earliest PubMed `ArticleDate`, otherwise the journal issue
  date. Records without either date are retained in the raw ledger but excluded
  from temporal windows.

## Retrieval, formula identity, and screening ledger

Only exact English aliases are queried in Title/Abstract. The broad aliases
`Rehmannia six`, `Four Agents decoction`, and `Xiaoyao powder` are not used.
The ledger retains every machine-retrieved record, including exclusions. An
automated proxy may use only records that have:

- no review/editorial/letter/news status;
- no retraction, correction, corrigendum, erratum, or protocol notice;
- an exact formula alias in the title; and
- no machine-detected modified/combined-formula marker in the title.

Records failing these rules remain available for blinded human review but cannot
contribute a machine proxy label. Manual review, not automatic triage, decides
whether a paper studies the exact formulation and whether a target is supported.

## Entity and pharmacology rules

Automated gene mentions are high-specificity screening cues. An exact HGNC
symbol match is necessary but not sufficient: known ambiguous short strings are
suppressed and other short symbols require local biological-target context.
Human full-text adjudication remains the reference standard. Current HGNC and
STP resources are not treated as historical snapshots. Current STP is a
non-temporal comparator and cannot expand the primary candidate universe.

ChEMBL activities are queried in the current database with `document_year <=
2022`, human-target and non-missing-pChEMBL filters. This is called
document-year constrained, not a historical ChEMBL snapshot. Missing herbs or
marker compounds are recorded as unavailable; they are never imputed as zero.

## Outcomes

The current manuscript reports retrieval counts, data-quality exclusions,
provenance, missingness, deterministic rankings, and software checks. It does
not report recall, AP, NDCG, superiority, calibration, or a biological-truth
claim. If a later gold-standard comparison is performed, it will use dual
independent screening, third-reviewer adjudication, PMID-clustered sensitivity,
multiple K values, and all-formula sensitivity reporting.
