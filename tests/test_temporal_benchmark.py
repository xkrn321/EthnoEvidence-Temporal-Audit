import math

from ethno_evidence.eval.temporal_benchmark import (
    assess_record_eligibility,
    average_precision_at_k,
    bootstrap_macro_ci,
    build_pubmed_query,
    gene_mentions_by_article,
    ndcg_at_k,
    recall_at_k,
    suggest_gene_evidence,
    triage_article,
)
from ethno_evidence.evidence.literature_mining import FORMULA_QUERIES


def test_query_is_time_bounded_and_alias_scoped():
    q = build_pubmed_query(["Sini Decoction", "Sini Tang"], 2023, 2026)
    assert '"Sini Decoction"[Title/Abstract]' in q
    assert '"2023/01/01"[Date - Publication]' in q
    assert '"2026/12/31"[Date - Publication]' in q
    assert '"Review"[Publication Type]' in q


def test_triage_excludes_review_and_prioritises_assays():
    review = {"publication_types": ["Review"], "title": "", "abstract": ""}
    assert triage_article(review)["status"] == "exclude"
    original = {
        "publication_types": ["Journal Article"],
        "title": "Mechanistic study",
        "abstract": "Western blot and qPCR were performed in mice.",
    }
    triage = triage_article(original)
    assert triage["status"] == "candidate_experimental"
    assert "protein_assay" in triage["reasons"]


def test_proxy_eligibility_excludes_status_and_formula_variants():
    aliases = ["Sijunzi Decoction", "Sijunzi Tang"]
    retracted = {"title": "Retracted: Sijunzi Decoction study",
                 "abstract": "", "publication_types": ["Journal Article"]}
    assert assess_record_eligibility(retracted, aliases)["status"] == "exclude"
    modified = {"title": "Modified Sijunzi Decoction improved a model",
                "abstract": "", "publication_types": ["Journal Article"]}
    assert assess_record_eligibility(modified, aliases)["status"] == "exclude"
    prefixed = {"title": "Huangqi Sijunzi Decoction improved a model",
                "abstract": "", "publication_types": ["Journal Article"]}
    assert assess_record_eligibility(prefixed, aliases)["status"] == "exclude"
    combined = {"title": "Sijunzi Decoction and Guipi Tang relieved a model",
                "abstract": "", "publication_types": ["Journal Article"]}
    assert assess_record_eligibility(combined, aliases)["status"] == "exclude"
    exact = {"title": "Effects of Sijunzi Decoction in a model",
             "abstract": "", "publication_types": ["Journal Article"]}
    assert assess_record_eligibility(exact, aliases)["status"] == "eligible"


def test_proxy_gene_mentions_require_specificity_and_context():
    records = {
        "1": {
            "title": "Metabolic assessment",
            "abstract": "TG fell in serum. AKT1 protein expression was reduced.",
            "publication_types": ["Journal Article"],
        }
    }
    out = gene_mentions_by_article(records, {"TG", "AKT1"},
                                   aliases=["Sijunzi Decoction"])
    assert out["1"]["genes"] == {"AKT1": 1}


def test_papera_aliases_remove_broad_legacy_terms():
    assert "Rehmannia six" not in FORMULA_QUERIES["六味地黄丸"]
    assert "Four Agents decoction" not in FORMULA_QUERIES["四物汤"]
    assert "Xiaoyao powder" not in FORMULA_QUERIES["逍遥散"]


def test_ranking_metrics_and_bootstrap_are_deterministic():
    ranked = ["A", "B", "C", "D"]
    positives = {"B", "D"}
    assert recall_at_k(ranked, positives, 2) == 0.5
    assert average_precision_at_k(ranked, positives, 4) == 0.5
    assert 0 < ndcg_at_k(ranked, {"B": 2, "D": 1}, 4) < 1
    assert math.isnan(recall_at_k(ranked, set(), 2))
    a = bootstrap_macro_ci([0.2, 0.4, 0.8], n_boot=200, seed=7)
    b = bootstrap_macro_ci([0.2, 0.4, 0.8], n_boot=200, seed=7)
    assert a == b
    assert a["n"] == 3


def test_gene_evidence_suggestion_includes_context_but_is_not_gold():
    record = {
        "title": "AKT1 experiment",
        "abstract": "AKT1 knockdown reversed the phenotype. TNF was discussed.",
    }
    out = suggest_gene_evidence(record, {"AKT1", "TNF"})
    assert out["AKT1"]["suggested_grade"] == 2
    assert "knockdown" in out["AKT1"]["context"]
    assert out["TNF"]["suggested_grade"] == "U"
