"""Frozen temporal-benchmark utilities for Paper A.

This module separates candidate retrieval from the human reference-standard
screen. Automated triage is never treated as a gold label. The post-cutoff
target set must be confirmed by two reviewers before inferential results are
reported.
"""
from __future__ import annotations

import math
import random
import re
from statistics import mean
from typing import Dict, Iterable, List, Mapping, Optional, Sequence

REVIEW_TYPES = {
    "Review", "Systematic Review", "Meta-Analysis", "Editorial", "Letter",
    "Comment", "News",
}

# Retraction notices, corrections, protocols and modified/combined formulae are
# preserved in the screening ledger but never enter an automated proxy label.
_STATUS_PATTERNS = re.compile(
    r"\b(retracted|retraction notice|corrigendum|erratum|correction to|"
    r"expression of concern|study protocol|protocol for)\b", re.I)
_VARIANT_PATTERNS = re.compile(
    r"\b(modified|modified formula|addition(?:s)?|combined with|combination of|"
    r"co[- ]?administration|herb pair|herb couple|plus)\b", re.I)
_ALLOWED_PREFIX_WORDS = {
    "a", "an", "the", "of", "for", "by", "with", "on", "in", "from",
    "using", "use", "effects", "effect", "mechanism", "study", "treatment",
    "therapy", "against", "after", "and", "or", "to",
}

# Exact HGNC matching alone is not enough for short uppercase abbreviations in
# abstracts. These strings have repeatedly occurred as non-gene concepts in
# this corpus and are not emitted by the automated screening aid.
_AMBIGUOUS_SYMBOLS = {
    "CS", "EC", "GC", "MAP", "NPS", "PC", "SDS", "TG", "TST",
}
_SHORT_GENE_CONTEXT = re.compile(
    r"\b(gene|genes|protein|mrna|transcript|expression|expressed|"
    r"western blot|immunoblot|immunofluorescen\w*|elisa|qpcr|rt[- ]?pcr|"
    r"knock(?:down|out)|silenc\w*|overexpress\w*|crispr|sirna|shrna)\b",
    re.I,
)

_EXPERIMENTAL_PATTERNS = {
    "protein_assay": re.compile(
        r"\b(western blot|immunoblot|immunohistochemistry|immunofluorescen\w*|"
        r"ELISA|flow cytometr\w*)\b", re.I),
    "transcript_assay": re.compile(
        r"\b(qPCR|RT[- ]?PCR|RNA[- ]?seq|transcriptom\w*)\b", re.I),
    "perturbation": re.compile(
        r"\b(knock(?:down|out)|silenc\w*|overexpress\w*|inhibitor|agonist|"
        r"antagonist|siRNA|shRNA|CRISPR|rescue experiment)\b", re.I),
    "binding_activity": re.compile(
        r"\b(binding assay|surface plasmon resonance|SPR|isothermal titration|"
        r"IC50|EC50|Ki|Kd)\b", re.I),
    "biological_model": re.compile(
        r"\b(mice|mouse|rats?|zebrafish|cell line|cells|in vivo|in vitro|"
        r"randomi[sz]ed|patients?|cohort)\b", re.I),
}

_COMPUTATIONAL_PATTERNS = re.compile(
    r"\b(network pharmacology|molecular docking|in silico|bioinformatics)\b",
    re.I,
)


def build_pubmed_query(aliases: Sequence[str], start_year: int,
                       end_year: int) -> str:
    """Build the protocol-locked PubMed Title/Abstract query."""
    if not aliases:
        raise ValueError("at least one formula alias is required")
    names = " OR ".join(f'"{a}"[Title/Abstract]' for a in aliases)
    concepts = (
        "mechanism*[Title/Abstract] OR pathway*[Title/Abstract] OR "
        "gene*[Title/Abstract] OR protein*[Title/Abstract] OR "
        "pharmacolog*[Title/Abstract]"
    )
    exclusions = " OR ".join(f'"{p}"[Publication Type]' for p in sorted(REVIEW_TYPES))
    return (
        f"({names}) AND ({concepts}) AND "
        f'("{start_year}/01/01"[Date - Publication] : '
        f'"{end_year}/12/31"[Date - Publication]) NOT ({exclusions})'
    )


def assess_record_eligibility(record: Mapping[str, object],
                              aliases: Sequence[str]) -> dict:
    """Classify a machine-retrieved record for the *automated proxy only*.

    This deliberately conservative rule does not substitute for reviewer
    screening. It removes obvious status notices and non-exact formula records
    from proxy metrics while retaining them in the ledger with an explanation.
    """
    publication_types = set(record.get("publication_types") or [])
    excluded = sorted(publication_types & REVIEW_TYPES)
    title = str(record.get("title") or "")
    abstract = str(record.get("abstract") or "")
    if excluded:
        return {"status": "exclude", "reasons": ["publication_type:" + x
                                                   for x in excluded]}
    if _STATUS_PATTERNS.search(title):
        return {"status": "exclude", "reasons": ["publication_status_notice"]}
    if _VARIANT_PATTERNS.search(title):
        return {"status": "exclude", "reasons": ["modified_or_combined_formula"]}

    # A title-level exact alias is necessary for the automated proxy. Records
    # where an alias occurs only in the abstract remain for manual review but
    # do not create machine labels. This prevents an acronym or a related
    # formula from silently becoming an "exact formula" association.
    title_lower = title.lower()
    matched_alias = None
    match_start = None
    for alias in aliases:
        match = re.search(rf"(?<![a-z0-9]){re.escape(alias.lower())}(?![a-z0-9])",
                          title_lower)
        if match:
            matched_alias = alias
            match_start = match.start()
            break
    if not matched_alias:
        if any(alias.lower() in abstract.lower() for alias in aliases):
            return {"status": "manual_review",
                    "reasons": ["formula_alias_not_in_title"]}
        return {"status": "exclude", "reasons": ["formula_alias_not_relocated"]}

    prefix_words = re.findall(r"[a-z]+", title_lower[:match_start])
    if prefix_words and prefix_words[-1] not in _ALLOWED_PREFIX_WORDS:
        return {"status": "exclude",
                "reasons": ["title_prefixed_formula_variant", f"alias:{matched_alias}"]}
    alias_end = match_start + len(matched_alias)
    title_tail = title_lower[alias_end:]
    if re.search(r"\b(?:and|with|plus)\s+(?:[a-z-]+\s+){0,4}"
                 r"(?:decoction|tang|wan|pill)\b", title_tail):
        return {"status": "exclude",
                "reasons": ["formula_combination_in_title", f"alias:{matched_alias}"]}
    if re.search(r"\b(?:with\s+)?combined\s+decoction\b", title_tail):
        return {"status": "exclude",
                "reasons": ["formula_combination_in_title", f"alias:{matched_alias}"]}
    return {"status": "eligible", "reasons": [f"exact_title_alias:{matched_alias}"]}


def triage_article(record: Mapping[str, object],
                   aliases: Optional[Sequence[str]] = None) -> dict:
    """Machine triage for reviewer prioritisation, not a gold label."""
    if aliases is not None:
        eligibility = assess_record_eligibility(record, aliases)
        if eligibility["status"] == "exclude":
            return eligibility
    publication_types = set(record.get("publication_types") or [])
    excluded = sorted(publication_types & REVIEW_TYPES)
    if excluded:
        return {"status": "exclude", "reasons": ["publication_type:" + x
                                                    for x in excluded]}
    text = f"{record.get('title', '')} {record.get('abstract', '')}"
    signals = sorted(name for name, pattern in _EXPERIMENTAL_PATTERNS.items()
                     if pattern.search(text))
    computational = bool(_COMPUTATIONAL_PATTERNS.search(text))
    if signals:
        return {
            "status": "candidate_experimental",
            "reasons": signals + (["computational_component"] if computational else []),
        }
    if computational:
        return {"status": "manual_review", "reasons": ["computational_only_possible"]}
    return {"status": "manual_review", "reasons": ["no_abstract_assay_signal"]}


def _strict_gene_mentions(text: str, dictionary: set) -> dict:
    """Return high-specificity automatic gene mentions for screening only."""
    from ..evidence.literature_mining import extract_genes
    initial = extract_genes(text, dictionary, min_count=1)
    kept = {}
    for gene, count in initial.items():
        if gene in _AMBIGUOUS_SYMBOLS:
            continue
        contexts = re.findall(rf"[^.?!]*\b{re.escape(gene)}\b[^.?!]*[.?!]?", text,
                              flags=re.I)
        if len(gene) <= 3 and not any(_SHORT_GENE_CONTEXT.search(x)
                                      for x in contexts):
            continue
        kept[gene] = count
    return kept


def gene_mentions_by_article(records: Mapping[str, Mapping[str, object]],
                             dictionary: set,
                             aliases: Optional[Sequence[str]] = None) -> Dict[str, dict]:
    """High-specificity automatic screening mentions with PMID provenance.

    These mentions are a triage aid, not a gene-level reference standard.
    """
    out: Dict[str, dict] = {}
    for pmid, record in sorted(records.items()):
        text = f"{record.get('title', '')} {record.get('abstract', '')}"
        mentions = _strict_gene_mentions(text, dictionary)
        eligibility = (assess_record_eligibility(record, aliases)
                       if aliases is not None else None)
        out[pmid] = {"genes": mentions,
                     "triage": triage_article(record),
                     "eligibility": eligibility}
    return out


def suggest_gene_evidence(record: Mapping[str, object],
                          genes: Iterable[str]) -> Dict[str, dict]:
    """Return review aids (context + tentative grade) for mentioned genes.

    Suggestions are deliberately conservative and must not be copied into the
    reference standard without human confirmation from the abstract/full text.
    """
    text = f"{record.get('title', '')}. {record.get('abstract', '')}"
    sentences = [s.strip() for s in re.split(r"(?<=[.!?])\s+", text) if s.strip()]
    out = {}
    for gene in sorted(set(genes)):
        contexts = [s for s in sentences
                    if re.search(rf"\b{re.escape(gene)}\b", s)]
        if not contexts:
            continue
        local = " ".join(contexts[:2])
        direct = any(_EXPERIMENTAL_PATTERNS[key].search(local)
                     for key in ("perturbation", "binding_activity"))
        measured = any(_EXPERIMENTAL_PATTERNS[key].search(local)
                       for key in ("protein_assay", "transcript_assay"))
        grade = 2 if direct else (1 if measured else "U")
        out[gene] = {"suggested_grade": grade, "context": local[:600]}
    return out


def recall_at_k(ranked: Sequence[str], positives: Iterable[str], k: int) -> float:
    positives = set(positives)
    if not positives:
        return math.nan
    return len(set(ranked[:k]) & positives) / len(positives)


def average_precision_at_k(ranked: Sequence[str], positives: Iterable[str],
                           k: int) -> float:
    positives = set(positives)
    if not positives:
        return math.nan
    hits = 0
    total = 0.0
    for i, gene in enumerate(ranked[:k], start=1):
        if gene in positives:
            hits += 1
            total += hits / i
    return total / min(len(positives), k)


def ndcg_at_k(ranked: Sequence[str], relevance: Mapping[str, int], k: int) -> float:
    """NDCG with integer relevance grades (0=unlabelled/nonrelevant)."""
    def _dcg(values: Sequence[int]) -> float:
        return sum((2 ** rel - 1) / math.log2(i + 2)
                   for i, rel in enumerate(values))

    observed = [int(relevance.get(g, 0)) for g in ranked[:k]]
    ideal = sorted((int(x) for x in relevance.values()), reverse=True)[:k]
    denominator = _dcg(ideal)
    return _dcg(observed) / denominator if denominator else math.nan


def bootstrap_macro_ci(values: Sequence[float], n_boot: int = 10_000,
                       seed: int = 20260828) -> dict:
    """Formula-level non-parametric bootstrap CI for a macro mean."""
    clean = [float(v) for v in values if not math.isnan(float(v))]
    if not clean:
        return {"n": 0, "mean": None, "ci95": [None, None]}
    rng = random.Random(seed)
    samples = sorted(mean(rng.choice(clean) for _ in clean)
                     for _ in range(n_boot))
    lo = samples[int(0.025 * (n_boot - 1))]
    hi = samples[int(0.975 * (n_boot - 1))]
    return {"n": len(clean), "mean": mean(clean), "ci95": [lo, hi],
            "n_boot": n_boot, "seed": seed}
