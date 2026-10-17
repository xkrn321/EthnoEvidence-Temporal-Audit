# EthnoEvidence Temporal Audit v2.0.2

This is the frozen, publication-focused v2.0.2 release accompanying the
methodology manuscript *Evidence-Gated Target Prioritization for Multicomponent
Pharmacology: A Temporally Audited Cheminformatics Workflow*. The v2.0.2 tag
remains unchanged; the current `main` branch includes follow-up figure and audit
clarifications described below.

It makes one bounded claim: evidence provenance, time eligibility, explicit
missingness, and declared claim ceilings can be represented as auditable inputs
to multicomponent target prioritization.  It does **not** estimate biological
target accuracy, establish an intact-formula mechanism, or substitute for
preparation detection, in-vivo exposure, or selective perturbation.

## Reproduce the frozen release

```bash
python -m pip install -e '.[dev]'
python scripts/verify_papera_v2_release.py
pytest -q
```

The offline verifier never calls live services. It checks all 23 locked
support-file checksums, the 2022/2023 time boundary, the complete 2023–2025
held-out screening ledger, formula-specific train/test PMID disjointness, and
exact reconstruction of all declared rankings and ablations from serialised
feature objects.

The frozen v2.0.2 tag contains 11 released tests. Current `main` additionally
tests figure denominators and available-k reporting; its expected output is 17
passed tests and a successful verifier report. No controlled-data tests are
shipped.

The current `main` branch labels Figure 3 overlap as available-k
(`k=min(10, N_gated, N_baseline)`) and Figure 4 turnover as full Top-10 only.
The measured-only comparison therefore uses 13 formulas; Liuwei Dihuang Wan
has only three measured-supported candidates and is excluded from full Top-10
replacement summaries. The audit table records its available `k` explicitly.

## Release scope

Included:

- a fixed 14-formula temporal audit and literature-screening resource;
- document-year-constrained human-target ChEMBL records, locked rankings,
  audits, figures, source code, and checksum manifest;
- machine-readable rules that preserve exclusions and avoid treating automated
  proxy eligibility as a human reference standard.

An exploratory structure-matched transfer audit is also provided in
[`transfer_audit/`](transfer_audit/). It contains the frozen protocol, derived
summary tables, task-level decision audits, and source-dependence diagnostics for
40 additional small-molecule tasks. The audit is explicitly exploratory and does
not estimate biological accuracy or comparator superiority. Licensed DrugCentral
raw tables and the local retrieval cache are not redistributed. The
`evidence_chain_trace.json` file provides a candidate-level provenance example
for the Huanglian Jiedu Decoction contrast and records the claim ceiling for
each evidence chain.

Excluded:

- practitioner interviews, formula variants from fieldwork, specimen/voucher
  information, traditional-knowledge records, and identifiable or restricted
  material;
- biological outcome labels and any claim of predictive performance.

See [PUBLIC_RELEASE_BOUNDARY.md](PUBLIC_RELEASE_BOUNDARY.md) for the scope and
limits of this public release.
