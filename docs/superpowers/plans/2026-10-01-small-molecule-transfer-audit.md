# Small Molecule Transfer Audit Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build and audit an independent, public-data small-molecule target-prioritisation case without changing the frozen 14-formula benchmark or claiming biological accuracy.

**Architecture:** A separate `transfer` package selects a deterministic drug cohort, normalises dated ChEMBL and PubMed evidence into provenance objects, and produces three fixed-universe rankings. A retrieval script stores raw responses and hashes; an offline verifier reconstructs all reported outputs. The manuscript is revised only after those outputs pass local checks.

**Tech Stack:** Python 3.9+, existing project dependencies, ChEMBL REST API, NCBI E-utilities, JSON/TSV, pytest, python-docx for the manuscript.

**Spec:** `docs/superpowers/specs/2026-10-01-small-molecule-transfer-audit-design.md`

## Global Constraints

- Keep `benchmarks/papera_temporal_v2/`, its release checksums, and the 14-formula reported results unchanged.
- Primary cutoff: 31 December 2022; 2023–2025 records are descriptive only.
- Cohort strata: ATC first letter C, J, N; priority C then J then N; at most six molecules per stratum; order by SHA-256 of `EthnoEvidence-transfer-2026-10-01|CHEMBL_ID`.
- ChEMBL molecule eligibility: `max_phase = 4`, `molecule_type = Small molecule`, `structure_type = MOL`, `first_approval` in 2010–2022, nonempty canonical SMILES, and self-parented identifier.
- Primary measured tier: human single-protein target, nonmissing pChEMBL, interpretable document and assay identifiers, `document_year <= 2022`, pChEMBL >= 5.0.
- Never turn an absent record, ambiguous identity, failed request, or machine-screening mention into a verified negative or biological gold label.
- Do not run new SwissTargetPrediction jobs or claim a historical ChEMBL snapshot from document-year filtering.
- Do not publish new data or push repository changes until author review and redistribution checks.

## Review Focus

1. **Salt/parent duplicates:** one parent appears once in the cohort and activity mapping; test in Task 1.
2. **Ambiguous names and short gene symbols:** retain manual-review status instead of automatic positive literature evidence; test in Task 2.
3. **Date gaps and current database changes:** missing document year is ineligible with a reason, and retrieval manifest records release/date/hash; test in Task 2.
4. **Empty candidate universe:** retain the drug in the denominator as not rankable, without a fabricated Top-10; test in Task 3.
5. **API interruption or pagination change:** no partial retrieval is labelled complete; test in Task 1.

---

### Task 1: Lock cohort selection and raw retrieval

**Files:**
- Create: `src/ethno_evidence/transfer/cohort.py`
- Create: `src/ethno_evidence/transfer/retrieve.py`
- Create: `benchmarks/small_molecule_transfer/protocol.json`
- Create: `scripts/build_small_molecule_transfer.py`
- Test: `tests/test_transfer_cohort.py`

**Interfaces:**
- Produces `select_cohort(molecules: Sequence[Mapping[str, object]]) -> list[dict]`, with each selected record containing `drug_id`, `stratum`, `parent_id`, `preferred_name`, `canonical_smiles`, and `first_approval`.
- Produces `fetch_all_molecules(session, *, page_size: int = 1000) -> tuple[list[dict], list[dict]]`; the second list is the page manifest. An incomplete page sequence raises rather than returning a partial cohort.
- Script CLI accepts `--retrieve` to fetch public data and `--offline` to rebuild from stored raw responses; both write under `benchmarks/small_molecule_transfer/` only.

- [ ] **Step 1: Write failing tests.** Assert all exact eligibility fields and the C/J/N priority; a self-parented duplicate appears once; six per stratum are chosen by the specified hash order; a short stratum is not topped up from another; interrupted pagination raises.
- [ ] **Step 2: Run `pytest -q tests/test_transfer_cohort.py` and confirm the tests fail for missing interfaces.**
- [ ] **Step 3: Implement the interfaces and commit the protocol before any full retrieval.** Raw responses, URL, retrieval UTC time, page metadata, and SHA-256 hashes go into a manifest; a failed page leaves status `incomplete`.
- [ ] **Step 4: Run `pytest -q tests/test_transfer_cohort.py` and require all tests to pass.**
- [ ] **Step 5: Commit only Task 1 code, tests, and the locked protocol.** Do not commit a live result at this step.

### Task 2: Normalise measured and literature evidence

**Files:**
- Create: `src/ethno_evidence/transfer/evidence.py`
- Create: `src/ethno_evidence/transfer/pubmed.py`
- Test: `tests/test_transfer_evidence.py`
- Modify: `scripts/build_small_molecule_transfer.py`

**Interfaces:**
- Consumes the selected records from Task 1.
- Produces `normalise_activity(row: Mapping[str, object], drug_id: str, cutoff_year: int = 2022) -> dict` with `target_id`, `document_id`, `assay_id`, `document_year`, `pchembl`, `eligibility`, and `reason`.
- Produces `screen_literature(record: Mapping[str, object], drug: Mapping[str, object], gene_dictionary: set[str], cutoff_year: int = 2022) -> dict` with PMID, canonical date, source status, target mentions, and manual-review reasons.
- Produces `build_evidence_objects(drug: Mapping[str, object], activities: Sequence[Mapping[str, object]], articles: Sequence[Mapping[str, object]]) -> list[dict]` in a documented JSON schema with source class, identity, date, eligibility, missingness reason, and claim ceiling.

- [ ] **Step 1: Write failing tests.** Include a valid human single-protein activity; nonhuman, non-single-protein, absent pChEMBL, missing year, and post-cutoff exclusions; ambiguous drug name and short gene symbol routed to manual review; API failure represented as missingness, not zero.
- [ ] **Step 2: Run `pytest -q tests/test_transfer_evidence.py` and confirm missing-interface failures.**
- [ ] **Step 3: Implement normalisation and a bounded PubMed retrieval with page/query/date metadata.** Reuse existing date and gene-screening utilities only where semantics match; save raw API responses and explicit exclusions.
- [ ] **Step 4: Run `pytest -q tests/test_transfer_evidence.py` and inspect a small offline fixture by hand.**
- [ ] **Step 5: Commit Task 2 code and tests.** No new manuscript claims.

### Task 3: Fixed-universe rankings and sensitivity audit

**Files:**
- Create: `src/ethno_evidence/transfer/ranking.py`
- Create: `benchmarks/small_molecule_transfer/RANKING_SPEC.md`
- Test: `tests/test_transfer_ranking.py`
- Modify: `scripts/build_small_molecule_transfer.py`

**Interfaces:**
- Consumes evidence objects from Task 2.
- Produces `make_target_features(drug_id: str, objects: Sequence[Mapping[str, object]]) -> dict[str, dict]` and `rank_transfer(features: Mapping[str, Mapping[str, object]], *, threshold: float = 5.0) -> dict[str, list[str]]` with keys `evidence_gated`, `measured_only`, and `ungated_union`.
- Produces `summarise_rank_changes(rankings: Mapping[str, Sequence[str]], k: int = 10) -> dict` with actual available-K denominators, overlap, replacements, and non-rankable status.

- [ ] **Step 1: Write failing tests.** Check one fixed union candidate universe, no formula cross-compound term, deterministic target-ID tie-break, thresholds 4.0/4.5/5.0/5.5/6.0, empty universe, and a measured-only list shorter than ten.
- [ ] **Step 2: Run `pytest -q tests/test_transfer_ranking.py` and confirm missing-interface failures.**
- [ ] **Step 3: Write and commit `RANKING_SPEC.md` with exact rank keys and missing-value behaviour before running the real cohort.** Confirm it does not redefine the original formula ranking.
- [ ] **Step 4: Implement the ranker and audit summary.** Keep all three outputs on the same fixed universe; in measured-only ranking, nonmeasured candidates follow measured candidates with explicit status rather than silently disappearing.
- [ ] **Step 5: Run `pytest -q tests/test_transfer_ranking.py` and require all tests to pass.**
- [ ] **Step 6: Commit Task 3 code and tests.**

### Task 4: Run and verify the independent cohort

**Files:**
- Create: `scripts/verify_small_molecule_transfer.py`
- Create: `benchmarks/small_molecule_transfer/RETRIEVAL_MANIFEST.json`
- Create: `benchmarks/small_molecule_transfer/cohort.json`
- Create: `benchmarks/small_molecule_transfer/evidence_objects.json`
- Create: `benchmarks/small_molecule_transfer/rankings.json`
- Create: `benchmarks/small_molecule_transfer/audit_summary.json`
- Test: `tests/test_transfer_release.py`

**Interfaces:**
- The verifier reads only stored files, recomputes hashes and selection, validates date eligibility, and exactly reconstructs rank lists and audit counts without network access.

- [ ] **Step 1: Write a failing verifier test** that detects a modified raw-response hash, a changed cohort member, a post-cutoff primary record, and an altered rank list.
- [ ] **Step 2: Run `pytest -q tests/test_transfer_release.py` and confirm the failures.**
- [ ] **Step 3: Implement the offline verifier and perform the frozen public retrieval.** Log incomplete API attempts; never use a partial data set as if complete. Record cold and cached elapsed time, platform, Python and dependency versions, and network-call counts.
- [ ] **Step 4: Run the new verifier, `pytest -q`, and the existing `python scripts/verify_papera_v2_release.py`.** Compare key outputs manually with raw records and record any failed transfer or empty-stratum result.
- [ ] **Step 5: Audit third-party redistribution terms before adding raw records to a public release.** If uncertain, keep local raw files and publish retrieval instructions plus hashes only.
- [ ] **Step 6: Commit verified Task 4 code and permitted outputs locally.** Do not push.

### Task 5: Update the manuscript only from verified results

**Files:**
- Modify: `../Frontiers_Research_Methods_Review/Manuscript_Research_Methods_Candidate.docx`
- Modify: `../Frontiers_Research_Methods_Review/Manuscript_Review.md`
- Create: `../Frontiers_Research_Methods_Review/Transfer_Claim_Ledger.md`

**Interfaces:**
- Consumes the verified Task 4 outputs and a human-reviewed claim ledger; produces an updated Word candidate, not a submission-approved manuscript.

- [ ] **Step 1: Create a claim ledger mapping every new number to a specific output field and retrieval manifest.** Mark weak coverage or failed transfer as such.
- [ ] **Step 2: Add a bounded exploratory transfer subsection and supplementary detail only if the design's decision rule is met.** Otherwise add the attempted transfer and its failure modes to Discussion; do not manufacture a positive validation story.
- [ ] **Step 3: Reconcile title, abstract, methods, figures/tables, discussion, and supplement with the claim ledger.** Keep biological accuracy and cross-disciplinary transfer explicitly untested.
- [ ] **Step 4: Render every Word page, visually inspect each page, check figure/table citations, and rerun numeric consistency checks.**
- [ ] **Step 5: Hand the revised candidate and remaining author-only declarations to the authors for scientific approval.** Do not submit or push without that review.
