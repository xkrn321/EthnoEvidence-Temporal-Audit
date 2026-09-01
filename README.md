# EthnoEvidence Temporal Audit v2.0.2

This is the frozen, publication-focused release accompanying the methodology
manuscript *Evidence-Gated Target Prioritization for Multicomponent
Pharmacology: A Temporally Audited Cheminformatics Workflow*.

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

Expected output is 11 passed tests and a successful verifier report. This
candidate intentionally retains only the tests that exercise the released
temporal methodology; no controlled-data tests are shipped.

## Release scope

Included:

- a 14-formula temporal calibration and literature-screening resource;
- document-year-constrained human-target ChEMBL records, locked rankings,
  audits, figures, source code, and checksum manifest;
- machine-readable rules that preserve exclusions and avoid treating automated
  proxy eligibility as a human reference standard.

Excluded:

- practitioner interviews, formula variants from fieldwork, specimen/voucher
  information, traditional-knowledge records, and identifiable or restricted
  material;
- biological outcome labels and any claim of predictive performance.

See [PUBLIC_RELEASE_BOUNDARY.md](PUBLIC_RELEASE_BOUNDARY.md) before making a
public archival release.
