# Small molecule transfer audit design

Status: design for author review. No cohort has been selected or analysed under this specification.

## Purpose and claim boundary

Test whether the provenance, temporal-eligibility, missingness, and claim-control components of EthnoEvidence can be applied to an independently assembled single-small-molecule target-prioritisation task. This is a transfer stress test, not a test of biological target accuracy, clinical utility, or cross-disciplinary generalisation. The existing 14-formula benchmark and its frozen release remain unchanged.

The scientific question is: when the same declared evidence classes and cutoff are applied to an independent drug cohort, can the workflow produce an auditable candidate list, and how much do the candidate lists change under alternative evidence rules? A null ranking change or a failed transfer is a reportable result.

## Design options considered

1. **Recommended: independent public-data transfer audit.** Assemble a deterministic cohort of approved single-molecule drugs from a frozen ChEMBL retrieval, add date-eligible measured and literature evidence, and run fixed-universe ranking comparisons. This directly tests reuse conditions but requires an adapter and careful identity checks.
2. **Documentation-only example.** Process one familiar drug to illustrate inputs and outputs. This is quick and useful for a tutorial but does not address the prior single-case criticism.
3. **Cross-disciplinary analogy.** Map evidence classes to a non-pharmacology decision problem. This would broaden the claim but would redefine the method and risk a superficial analogy; it is excluded from this round.

## Cohort and data policy

The cohort will be selected before looking at ranking outcomes. Retrieve all ChEMBL molecules meeting `max_phase = 4`, `molecule_type = Small molecule`, `structure_type = MOL`, `first_approval` from 2010 through 2022 inclusive, a nonempty canonical SMILES, and a self-parented molecule identifier. Use the first ATC letter to assign three prespecified strata: C (cardiovascular), J (systemic anti-infective), and N (nervous system). A molecule with multiple eligible strata is assigned C, then J, then N. Within each stratum, order identifiers by SHA-256 of `EthnoEvidence-transfer-2026-10-01|CHEMBL_ID` and take the first six. If a stratum contains fewer than six eligible molecules, include all and report the shortfall. Do not replace a selected drug because it has few activities, little literature, or an inconvenient result. Record all initial counts, exclusions, API URLs, retrieval date, ChEMBL release information when available, and SHA-256 hashes of raw responses.

The analytic unit is one selected drug. The fixed candidate universe for that drug is the union of eligible human protein targets appearing in either the measured or literature layer before the cutoff. A drug with an empty universe remains in the cohort and is reported as not rankable; it is not silently replaced. The primary temporal cutoff is 31 December 2022, matching the existing benchmark, and 2023–2025 records are reserved for a descriptive audit only. Current database contents filtered by document year are not described as a reconstructed historical database snapshot.

## Evidence layers and identity

- **Measured:** ChEMBL activity rows linked to the selected parent molecule, Homo sapiens single-protein targets, nonmissing pChEMBL, interpretable assay and document identifiers, and `document_year <= 2022`. Preserve all retrieved rows and exclusion reasons. A pChEMBL threshold of 5.0 defines the primary qualifying tier; thresholds 4.0, 4.5, 5.5, and 6.0 are sensitivity checks, not tuning choices.
- **Literature:** PubMed records retrieved with the exact preferred drug name and explicitly listed unambiguous synonyms, restricted by publication date. HGNC mentions are screening aids, not experimentally established drug–target relations. Preserve PMID, title, abstract availability, publication/date fields, query, and exclusion reasons. Ambiguous drug names or gene symbols are marked for manual review; they do not automatically become positive evidence.
- **Prediction:** no new SwissTargetPrediction jobs or other current prediction outputs are used in the primary historical comparison. If an extant prediction layer is explored later, it must remain explicitly non-temporal and cannot add primary candidates.
- **Missingness:** unresolved identity, inaccessible records, absent assays, absent abstracts, API failures, and unmapped targets are distinct statuses, never biological negatives.

## Adapter and ranking behaviour

Create an isolated small-molecule adapter that converts source rows into versioned evidence objects with `drug_id`, `target_id`, `source_class`, `source_record_id`, `source_date`, `eligibility`, `missingness_reason`, and `claim_ceiling`. Do not change the frozen formula benchmark, its ranking files, or its figure data. Reuse the existing temporal and provenance checks where their semantics match; document every rule that is not portable.

For each rankable drug, generate three deterministic orders over exactly the same candidate universe: evidence-gated, measured-only, and ungated union. The single-compound adapter removes the formula-specific cross-compound coverage term because it is constant, and records that deviation rather than claiming the original formula ranking is transferred unchanged. Rank keys, tie-breaks, and missing-value treatment must be written to a machine-readable specification before running the cohort. Do not optimise weights or thresholds against the observed Top-10 changes.

## Outputs and success criteria

Publish a cohort flow table, per-drug source coverage and missingness ledger, target-level provenance objects, rank lists, Top-10 overlap/replacement distributions, threshold sensitivity, and failure reasons. Report the number of rankable drugs and the distribution across strata. Use descriptive statistics only; these drugs are not independent samples from all pharmacology and there is no target-truth reference standard. A successful *feasibility* result means that at least one drug can be processed end to end with an inspectable evidence trail and deterministic reconstruction. A stronger transfer claim requires the method to run without case-specific manual correction on most selected drugs; the denominator and every exception must be shown. Neither outcome establishes biological superiority.

## Verification and release boundaries

Tests must cover deterministic selection, duplicate/parent handling, date eligibility, target identity, missingness, fixed-universe comparators, and exact rerun of rankings from serialized feature objects. Run the existing release verifier to confirm no regression. Record elapsed time and resource use for cold and cached runs, with operating system, Python version, dependency versions, and whether network calls were made. Check public-data redistribution terms before adding raw third-party records to a release. The manuscript may describe the case only after the data, outputs, and tests are locally audited; the public repository is updated only after author review.

## Manuscript decision rule

If the adapter produces auditable outputs with meaningful coverage, add a clearly labelled exploratory transfer subsection and supplement, and revise the abstract without claiming external biological validation. If coverage is poor or major rules require redesign, report the attempted transfer and failure modes in Discussion or supplementary methods rather than adding a positive-looking second benchmark. Maintain the original 14-formula results and declarations unchanged except for verified corrections. Human authors must approve all new scientific interpretations and submission declarations.
