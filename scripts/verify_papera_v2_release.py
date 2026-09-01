"""Offline integrity audit for the frozen temporal resource v2.

This verifier deliberately does not call PubMed, ChEMBL, or
SwissTargetPrediction.  It checks the local release inputs that support the
accompanying manuscript: the checksum manifest, temporal split, screening rule,
and all locked rankings.  It is not a performance evaluation and cannot turn
the automated proxy into a human reference standard.
"""
from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import random
import sys
from pathlib import Path
from typing import Any, Callable


ROOT = Path(__file__).resolve().parents[1]
BENCH = ROOT / "benchmarks" / "papera_temporal_v2"
sys.path.insert(0, str((ROOT / "src").resolve()))

from ethno_evidence.eval.temporal_benchmark import (  # noqa: E402
    assess_record_eligibility,
    triage_article,
)


CORE_FILES = (
    "protocol.json",
    "PROTOCOL.md",
    "RANKING_SPEC.md",
    "DEVIATIONS.md",
    "candidate_corpus.json",
    "manual_screening.tsv",
    "cutoff_chembl_2022.json",
    "locked_rankings_and_proxy_audit.json",
    "retrieval_manifest.json",
)
RANKING_NAMES = {
    "evidence_gated",
    "pre_cutoff_literature_frequency",
    "cutoff_chembl_measured_only",
    "stp_only_non_temporal_comparator",
    "ungated_union",
    "seeded_random",
    "without_measured_tier",
    "without_pre_cutoff_literature",
    "without_missingness_penalty",
}


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def _read_json(name: str) -> dict[str, Any]:
    return json.loads((BENCH / name).read_text(encoding="utf-8"))


def _ranker() -> Any:
    path = ROOT / "scripts" / "build_temporal_rankings.py"
    spec = importlib.util.spec_from_file_location("papera_locked_ranker", path)
    _require(spec is not None and spec.loader is not None,
             "could not load locked ranking builder")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _verify_checksums() -> int:
    lock = _read_json("LOCK.sha256.json")
    expected_paths = list(CORE_FILES) + [
        str(path.relative_to(BENCH))
        for path in sorted((BENCH / "raw_pubmed").glob("*_all.json"))
    ]
    entries = lock.get("files", [])
    _require(lock.get("release_id") == "papera_temporal_v2",
             "checksum manifest has the wrong release ID")
    _require(lock.get("checksum_algorithm") == "sha256",
             "checksum manifest has the wrong algorithm")
    _require(lock.get("file_count") == len(expected_paths),
             "checksum manifest file count differs from release inputs")
    _require([entry.get("path") for entry in entries] == expected_paths,
             "checksum manifest does not cover exactly the release inputs")
    for entry in entries:
        path = BENCH / entry["path"]
        _require(path.is_file(), f"checksum input is missing: {entry['path']}")
        _require(_sha256(path) == entry["sha256"],
                 f"checksum mismatch: {entry['path']}")
    return len(entries)


def _verify_temporal_split() -> tuple[int, int, int]:
    protocol = _read_json("protocol.json")
    manifest = _read_json("retrieval_manifest.json")
    corpus = _read_json("candidate_corpus.json")
    locked = _read_json("locked_rankings_and_proxy_audit.json")
    _require(protocol["protocol_id"] == manifest["protocol_id"] ==
             corpus["protocol_id"] == locked["protocol_id"] ==
             "papera_temporal_v2", "protocol IDs are inconsistent")
    _require(protocol["evaluation_start_year"] == 2023 and
             protocol["evaluation_end_year"] == 2025,
             "primary endpoint is not the complete 2023-2025 window")
    _require("2026" not in manifest["primary_evaluation_window"],
             "primary evaluation window includes incomplete 2026")
    _require(manifest["canonical_date_policy"] ==
             "first_article_date_else_journal_issue_date",
             "canonical publication-date policy is inconsistent")
    _require(manifest["train_test_disjoint_within_formula"] is True,
             "manifest does not attest a disjoint temporal split")

    formula_names = set(protocol["formulas"])
    _require(set(corpus["formulas"]) == formula_names == set(locked["formulas"]),
             "formula registries differ across locked files")
    rows = list(csv.DictReader((BENCH / "manual_screening.tsv").open(
        encoding="utf-8"), delimiter="\t"))
    row_index = {(row["formula"], row["pmid"]): row for row in rows}
    _require(len(row_index) == len(rows), "screening ledger repeats a formula-PMID row")

    train_total = test_total = proxy_total = 0
    for formula in protocol["formulas"]:
        entry = corpus["formulas"][formula]
        raw_path = BENCH / "raw_pubmed" / f"{formula}_all.json"
        raw = json.loads(raw_path.read_text(encoding="utf-8"))
        _require(raw["formula"] == formula, f"raw PubMed formula mismatch: {formula}")
        _require("2026/" not in raw["query"],
                 f"raw retrieval has an incomplete 2026 window: {formula}")
        records = raw["records"]
        train = set(entry["train"]["pmids"])
        test = set(entry["test"]["pmids"])
        _require(len(train) == len(entry["train"]["pmids"]),
                 f"duplicate train PMID in {formula}")
        _require(len(test) == len(entry["test"]["pmids"]),
                 f"duplicate test PMID in {formula}")
        _require(not train & test, f"train-test PMID overlap in {formula}")
        _require(train | test <= set(records),
                 f"temporal PMID absent from raw ledger in {formula}")

        for pmid in train:
            record = records[pmid]
            _require(record["publication_date_policy"] ==
                     "first_article_date_else_journal_issue_date",
                     f"wrong date policy for train PMID {pmid}")
            _require(int(record["canonical_publication_year"]) <= 2022,
                     f"post-cutoff PMID in train split: {formula}/{pmid}")
        test_rows = {(formula, pmid) for pmid in test}
        _require(test_rows <= set(row_index),
                 f"screening ledger lacks a held-out PMID for {formula}")
        for pmid in test:
            record = records[pmid]
            row = row_index[(formula, pmid)]
            year = int(record["canonical_publication_year"])
            _require(2023 <= year <= 2025,
                     f"held-out PMID outside complete primary window: {formula}/{pmid}")
            _require(row["canonical_publication_date"] ==
                     record["canonical_publication_date"],
                     f"canonical date differs between raw and screen: {formula}/{pmid}")
            _require(int(row["year"]) == year,
                     f"canonical year differs between raw and screen: {formula}/{pmid}")

        expected_proxy = {
            pmid for pmid in test
            if (assess_record_eligibility(records[pmid], entry["aliases"])["status"] ==
                "eligible" and
                triage_article(records[pmid])["status"] == "candidate_experimental")
        }
        actual_proxy = set(entry["test"]["proxy_eligible_pmids"])
        _require(actual_proxy == expected_proxy,
                 f"proxy eligibility does not match the strict rule: {formula}")
        for pmid in test:
            row = row_index[(formula, pmid)]
            rule = assess_record_eligibility(records[pmid], entry["aliases"])
            triage = triage_article(records[pmid])
            _require(row["proxy_eligibility"] == rule["status"],
                     f"screening status differs from strict rule: {formula}/{pmid}")
            _require(row["auto_triage"] == triage["status"],
                     f"screening triage differs from the strict rule: {formula}/{pmid}")
            _require((pmid in actual_proxy) ==
                     (rule["status"] == "eligible" and
                      triage["status"] == "candidate_experimental"),
                     f"inconsistent proxy PMID membership: {formula}/{pmid}")
        _require(set(locked["formulas"][formula]["automated_proxy_eligible_pmids"]) ==
                 actual_proxy, f"locked proxy PMID list differs from corpus: {formula}")
        train_total += len(train)
        test_total += len(test)
        proxy_total += len(actual_proxy)

    _require(set(row_index) == {
        (formula, pmid)
        for formula, entry in corpus["formulas"].items()
        for pmid in entry["test"]["pmids"]
    }, "screening ledger has rows outside the held-out corpus")
    _require(train_total == manifest["train_article_count"],
             "train article count disagrees with manifest")
    _require(test_total == manifest["test_article_count"] == len(rows),
             "test article count disagrees with manifest or screening ledger")
    _require(proxy_total == manifest["test_proxy_eligible_article_count"],
             "proxy article count disagrees with manifest")
    return train_total, test_total, proxy_total


def _as_features(raw: dict[str, dict[str, Any]]) -> dict[str, dict[str, Any]]:
    return {
        gene: {
            **feature,
            "strong_compounds": set(feature["strong_compounds"]),
            "weak_compounds": set(feature["weak_compounds"]),
            "stp_compounds": set(feature["stp_compounds"]),
            "literature_pmids": set(feature["literature_pmids"]),
        }
        for gene, feature in raw.items()
    }


def _verify_rankings() -> int:
    protocol = _read_json("protocol.json")
    locked = _read_json("locked_rankings_and_proxy_audit.json")
    ranker = _ranker()
    for formula in protocol["formulas"]:
        result = locked["formulas"][formula]
        features = _as_features(result["features"])
        primary = {
            gene: feature for gene, feature in features.items()
            if (feature["strong_compounds"] or feature["weak_compounds"] or
                feature["literature_pmids"])
        }
        _require(len(primary) == result["candidate_count"],
                 f"candidate count is not the primary candidate universe: {formula}")
        rankings = result["rankings"]
        _require(set(rankings) == RANKING_NAMES,
                 f"ranking set differs from the locked specification: {formula}")

        max_lit = max((len(f["literature_pmids"]) for f in primary.values()), default=0)
        max_cov = max((len(f["strong_compounds"] | f["weak_compounds"])
                       for f in primary.values()), default=0)
        seed = (int(hashlib.sha256(formula.encode()).hexdigest()[:8], 16) +
                int(protocol["bootstrap"]["seed"]))
        seeded = sorted(primary)
        random.Random(seed).shuffle(seeded)
        expected: dict[str, list[str]] = {
            "evidence_gated": ranker._sorted(primary, include_stp=False),
            "pre_cutoff_literature_frequency": sorted(
                (gene for gene, f in primary.items() if f["literature_pmids"]),
                key=lambda gene: (-len(primary[gene]["literature_pmids"]), gene)),
            "cutoff_chembl_measured_only": ranker._measured_only(primary),
            "stp_only_non_temporal_comparator": sorted(
                primary,
                key=lambda gene: (-int(bool(primary[gene]["stp_compounds"])),
                                  -primary[gene]["best_stp_probability"],
                                  -len(primary[gene]["stp_compounds"]), gene)),
            "ungated_union": sorted(
                primary,
                key=lambda gene: (-ranker._weighted_score(primary[gene], max_lit,
                                                           max_cov), gene)),
            "seeded_random": seeded,
            "without_measured_tier": ranker._sorted(primary, measured=False),
            "without_pre_cutoff_literature": ranker._sorted(primary, literature=False),
            "without_missingness_penalty": sorted(
                primary,
                key=lambda gene: (-ranker._weighted_score(primary[gene], max_lit,
                                                           max_cov, True), gene)),
        }
        for name, expected_order in expected.items():
            actual_order = rankings[name]
            _require(actual_order == expected_order,
                     f"locked ranking diverges from its declared algorithm: {formula}/{name}")
            _require(len(actual_order) == len(set(actual_order)),
                     f"duplicate gene within a ranking: {formula}/{name}")
            _require(set(actual_order) <= set(primary),
                     f"ranking contains a non-primary candidate: {formula}/{name}")
        _require(set(rankings["stp_only_non_temporal_comparator"]) == set(primary),
                 f"STP comparator changed the primary candidate universe: {formula}")
    return len(protocol["formulas"])


def _verify_runtime_dependency_declarations() -> None:
    pyproject = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    _require('"beautifulsoup4>=4.12"' in pyproject,
             "beautifulsoup4 is imported but not declared as a runtime dependency")
    _require('"matplotlib>=3.7"' in pyproject,
             "matplotlib is required for manuscript figures but not declared")


def main() -> int:
    try:
        _verify_runtime_dependency_declarations()
        checksum_count = _verify_checksums()
        train_total, test_total, proxy_total = _verify_temporal_split()
        formula_count = _verify_rankings()
    except (AssertionError, KeyError, TypeError, ValueError) as exc:
        print(f"[FAIL] temporal resource v2 release audit: {exc}", file=sys.stderr)
        return 1
    print("[OK] temporal resource v2 release audit")
    print(f"  checksums: {checksum_count} files")
    print(f"  temporal corpus: {formula_count} formulas; {train_total} pre-cutoff, "
          f"{test_total} held-out (2023-2025), {proxy_total} proxy-eligible")
    print("  rankings: all declared baselines and ablations reproduce exactly")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
