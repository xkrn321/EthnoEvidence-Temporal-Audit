"""Build locked temporal rankings and an explicitly non-gold pilot audit."""
from __future__ import annotations

import hashlib
import json
import math
import random
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str((ROOT / "src").resolve()))

from ethno_evidence.data.pubchem_client import PubChemClient  # noqa: E402
from ethno_evidence.data.swiss_client import SwissTargetPredictionClient  # noqa: E402
from ethno_evidence.eval.temporal_benchmark import (  # noqa: E402
    average_precision_at_k,
    bootstrap_macro_ci,
    recall_at_k,
)
from ethno_evidence.evidence.literature_mining import load_gene_dictionary  # noqa: E402
from ethno_evidence.evidence.reproduce import normalize  # noqa: E402

BENCH = ROOT / "benchmarks" / "papera_temporal_v2"


def _new_feature() -> dict:
    return {
        "strong_compounds": set(), "weak_compounds": set(),
        "stp_compounds": set(), "literature_pmids": set(),
        "best_pchembl": 0.0, "best_stp_probability": 0.0,
    }


def _rank_key(item: tuple[str, dict], measured: bool = True,
              literature: bool = True, include_stp: bool = False) -> tuple:
    gene, f = item
    strong = bool(f["strong_compounds"]) if measured else False
    weak = bool(f["weak_compounds"]) if measured else False
    stp = bool(f["stp_compounds"]) if include_stp else False
    lit = bool(f["literature_pmids"]) if literature else False
    tier = 4 if strong else 3 if weak else 2 if stp else 1 if lit else 0
    sources = sum((strong or weak, stp, lit))
    compounds = set()
    if measured:
        compounds |= f["strong_compounds"] | f["weak_compounds"]
    if include_stp:
        compounds |= f["stp_compounds"]
    return (tier, sources, len(compounds), f["best_pchembl"] if measured else 0,
            f["best_stp_probability"] if include_stp else 0,
            len(f["literature_pmids"]) if literature else 0)


def _weighted_score(f: dict, max_lit: int, max_coverage: int,
                    observed_only: bool = False,
                    include_stp: bool = False) -> float:
    measured = max(0.0, min(1.0, (f["best_pchembl"] - 4.0) / 4.0))
    lit = (math.log1p(len(f["literature_pmids"])) / math.log1p(max_lit)
           if max_lit else 0.0)
    coverage_n = len(f["strong_compounds"] | f["weak_compounds"] |
                     (f["stp_compounds"] if include_stp else set()))
    coverage = coverage_n / max_coverage if max_coverage else 0.0
    values = [measured, lit, coverage]
    present_flags = [
        bool(f["strong_compounds"] or f["weak_compounds"]),
        bool(f["literature_pmids"]),
        bool(coverage_n),
    ]
    if include_stp:
        values.insert(1, max(0.0, min(1.0, f["best_stp_probability"])))
        present_flags.insert(1, bool(f["stp_compounds"]))
    if observed_only:
        present = [value for value, yes in zip(values, present_flags) if yes]
        return sum(present) / len(present) if present else 0.0
    return sum(values) / len(values)


def _sorted(features: dict, measured: bool = True,
            literature: bool = True, include_stp: bool = False) -> list[str]:
    return [gene for gene, _ in sorted(
        features.items(), key=lambda item: (
            tuple(-x for x in _rank_key(item, measured, literature, include_stp)), item[0]))]


def _measured_only(features: dict) -> list[str]:
    """Pure ChEMBL-only baseline; STP/literature never break ties."""
    eligible = {
        gene: f for gene, f in features.items()
        if f["strong_compounds"] or f["weak_compounds"]
    }
    return sorted(eligible, key=lambda gene: (
        -int(bool(eligible[gene]["strong_compounds"])),
        -eligible[gene]["best_pchembl"],
        -len(eligible[gene]["strong_compounds"] | eligible[gene]["weak_compounds"]),
        gene,
    ))


def main() -> int:
    protocol = json.loads((BENCH / "protocol.json").read_text(encoding="utf-8"))
    chembl = json.loads((BENCH / "cutoff_chembl_2022.json").read_text(encoding="utf-8"))
    corpus = json.loads((BENCH / "candidate_corpus.json").read_text(encoding="utf-8"))
    approved = load_gene_dictionary()
    compounds = sorted(chembl["compounds"])

    pubchem = PubChemClient()
    name_smiles = {name: smiles for name in compounds
                   if (smiles := pubchem.get_smiles(name))}
    stp = SwissTargetPredictionClient().predict_many(name_smiles, workers=6)

    target_id_to_gene = {}
    for result in stp.values():
        for target in result.get("targets", []):
            gene = (target.get("common_name") or "").strip()
            target_id = target.get("chembl_id")
            if target_id and gene in approved:
                target_id_to_gene[target_id] = gene

    outputs = {"protocol_id": protocol["protocol_id"],
               "post_cutoff_gold_labels_present": False, "formulas": {}}
    pilot_metrics = defaultdict(list)
    for formula in protocol["formulas"]:
        features = defaultdict(_new_feature)
        formula_compounds = sorted({c for cs in chembl["formula_compounds"][formula].values()
                                    for c in cs})
        for compound in formula_compounds:
            measured = chembl["compounds"].get(compound, {})
            for tier in ("strong", "weak"):
                for target_id, row in measured.get(tier, {}).items():
                    gene = (normalize(row.get("target_pref_name") or "") or
                            target_id_to_gene.get(target_id))
                    if not gene or gene not in approved:
                        continue
                    features[gene][f"{tier}_compounds"].add(compound)
                    features[gene]["best_pchembl"] = max(
                        features[gene]["best_pchembl"], row["best_pchembl"])
            for target in stp.get(compound, {}).get("targets", []):
                gene = (target.get("common_name") or "").strip()
                if gene not in approved:
                    continue
                features[gene]["stp_compounds"].add(compound)
                features[gene]["best_stp_probability"] = max(
                    features[gene]["best_stp_probability"],
                    float(target.get("probability") or 0.0))
        train_mentions = corpus["formulas"][formula]["train"]["mentions"]
        for pmid, item in train_mentions.items():
            for gene in item["genes"]:
                features[gene]["literature_pmids"].add(pmid)

        # Current STP remains available as a comparator only. It cannot expand
        # the primary candidate universe or influence the evidence-gated rank.
        primary_features = {
            gene: f for gene, f in features.items()
            if f["strong_compounds"] or f["weak_compounds"] or f["literature_pmids"]
        }
        max_lit = max((len(x["literature_pmids"])
                       for x in primary_features.values()), default=0)
        max_cov = max((len(x["strong_compounds"] | x["weak_compounds"])
                       for x in primary_features.values()), default=0)
        gated = _sorted(primary_features, include_stp=False)
        measured_only = _measured_only(primary_features)
        stp_only = sorted(
            primary_features,
            key=lambda g: (-int(bool(primary_features[g]["stp_compounds"])),
                           -primary_features[g]["best_stp_probability"],
                           -len(primary_features[g]["stp_compounds"]), g))
        literature_only = sorted(
            (g for g, f in primary_features.items() if f["literature_pmids"]),
            key=lambda g: (-len(primary_features[g]["literature_pmids"]), g))
        ungated = sorted(primary_features, key=lambda g: (
            -_weighted_score(primary_features[g], max_lit, max_cov), g))
        no_missing_penalty = sorted(primary_features, key=lambda g: (
            -_weighted_score(primary_features[g], max_lit, max_cov, True), g))
        seed = (int(hashlib.sha256(formula.encode()).hexdigest()[:8], 16) +
                int(protocol["bootstrap"]["seed"]))
        seeded_random = sorted(primary_features)
        random.Random(seed).shuffle(seeded_random)
        rankings = {
            "evidence_gated": gated,
            "pre_cutoff_literature_frequency": literature_only,
            "cutoff_chembl_measured_only": measured_only,
            "stp_only_non_temporal_comparator": stp_only,
            "ungated_union": ungated,
            "seeded_random": seeded_random,
            "without_measured_tier": _sorted(primary_features, measured=False),
            "without_pre_cutoff_literature": _sorted(primary_features, literature=False),
            "without_missingness_penalty": no_missing_penalty,
        }

        # Automated post-cutoff mentions are an exploratory sensitivity audit,
        # never the locked gold standard.
        test_mentions = corpus["formulas"][formula]["test"]["mentions"]
        proxy_positives = {
            gene for item in test_mentions.values()
            if item["triage"]["status"] == "candidate_experimental"
            and item["eligibility"]["status"] == "eligible"
            for gene in item["genes"]
        }
        proxy = {}
        for method, ranking in rankings.items():
            rec = recall_at_k(ranking, proxy_positives, 20)
            ap = average_precision_at_k(ranking, proxy_positives, 20)
            # A formula with no eligible proxy gene has no denominator. Store
            # an explicit null in the JSON ledger rather than a non-standard
            # NaN, and do not treat it as a zero-performance result.
            proxy[method] = {
                "recall_at_20": None if math.isnan(rec) else rec,
                "ap_at_20": None if math.isnan(ap) else ap,
            }
            if not math.isnan(rec):
                pilot_metrics[(method, "recall_at_20")].append(rec)
            if not math.isnan(ap):
                pilot_metrics[(method, "ap_at_20")].append(ap)
        serial_features = {
            gene: {**f,
                   "strong_compounds": sorted(f["strong_compounds"]),
                   "weak_compounds": sorted(f["weak_compounds"]),
                   "stp_compounds": sorted(f["stp_compounds"]),
                   "literature_pmids": sorted(f["literature_pmids"])}
            for gene, f in primary_features.items()
        }
        outputs["formulas"][formula] = {
            "candidate_count": len(primary_features),
            "candidate_universe_policy": (
                "pre-cutoff ChEMBL and pre-cutoff literature only; current STP "
                "cannot add a primary candidate"),
            "formula_compound_count": len(formula_compounds),
            "features": serial_features,
            "rankings": rankings,
            "automated_proxy_eligible_pmids": corpus["formulas"][formula]["test"][
                "proxy_eligible_pmids"],
            "automated_proxy_positive_count": len(proxy_positives),
            "automated_proxy_metrics_not_gold": proxy,
        }

    summary = {}
    for (method, metric), values in pilot_metrics.items():
        summary.setdefault(method, {})[metric] = bootstrap_macro_ci(
            values,
            n_boot=int(protocol["bootstrap"]["replicates"]),
            seed=int(protocol["bootstrap"]["seed"]),
        )
    outputs["automated_proxy_summary_not_gold"] = summary
    path = BENCH / "locked_rankings_and_proxy_audit.json"
    path.write_text(json.dumps(outputs, ensure_ascii=False, indent=2,
                               allow_nan=False), encoding="utf-8")
    print(f"[OK] locked rankings -> {path}")
    for method, metrics in summary.items():
        r = metrics["recall_at_20"]
        print(f"{method:40s} proxy Recall@20={r['mean']:.3f} "
              f"[{r['ci95'][0]:.3f}, {r['ci95'][1]:.3f}]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
