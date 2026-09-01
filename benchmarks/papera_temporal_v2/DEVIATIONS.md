# Version change and deviations from v1

## 2026-08-29 integrity amendment

The v1 protocol claimed a 2023--2026 held-out window but relied on separate
PubMed ESearch windows and a parser that preferred `ArticleDate`. Thirteen
PMIDs therefore appeared in both train and test artefacts. The 2026 window was
also incomplete at the time of analysis. v2 supersedes v1 for all manuscript
work and makes the following pre-submission changes:

1. Retrieve one broad 1900--2025 range per formula, then assign windows with
   one canonical first-publication-date policy and enforce within-formula PMID
   disjointness.
2. Use complete calendar years 2023--2025 as the primary held-out screening
   window; 2026 is not in a primary endpoint.
3. Remove status notices, protocol notices, modified/combined formula records,
   and title-nonexact records from any automated proxy while retaining them in
   the screening ledger.
4. Replace raw HGNC-string matching in the automated screening aid with a
   high-specificity short-symbol/context rule.
5. Prevent current STP output from expanding the primary candidate universe or
   altering the primary evidence-gated rank.
6. Correct the pure measured-only comparator and delete redundant copies of the
   ungated ranking from the ablation set.

The amendment is not a post-hoc performance optimization: no gold labels are
available and no predictive-performance result is reported in the accompanying
manuscript.
