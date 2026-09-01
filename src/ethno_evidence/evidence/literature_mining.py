"""PubMed literature mining with temporal split (prospective validation).

Closes the loop between the static hand-written formula benchmarks and a
dynamic evidence base: ``lit_genes`` / ``lit_pathways`` are mined from
PubMed abstracts (E-utilities esearch + efetch), every gene appears
verbatim in a fetched abstract with PMID provenance, and a **temporal
split** (train <= cutoff / test > cutoff) is recorded so that any
prospective validation cannot be circular (predictions are judged
against genes reported only *after* the training window).

Honesty guards (iron rule 2 — no fabricated references):
- the HUGO dictionary is the default extraction source; symbols absent
  from an abstract are never included;
- the optional LLM-assisted step is *whitelist-filtered*: an LLM-suggested
  gene is kept only if it appears verbatim in the abstract text AND
  matches HUGO symbol syntax — fabrication is structurally impossible;
- every record receives an explicit first-publication-date policy; records
  without a parseable date are retained for audit but excluded from temporal
  validation windows.

Wording is restricted to "consistent with / reported in"; nothing here
claims mechanism confirmation.
"""
from __future__ import annotations

import json
import csv
import re
import time
import xml.etree.ElementTree as ET
from datetime import date
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import requests

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# 经典方 → 英文名（用于 PubMed 检索；无可靠英文名的方剂自动跳过）
FORMULA_QUERIES: Dict[str, List[str]] = {
    "逍遥散": ["Xiaoyao San", "Xiaoyaosan"],
    "四逆汤": ["Sini Decoction", "Sini Tang"],
    "六味地黄丸": ["Liuwei Dihuang Wan", "Liuwei Dihuang pill"],
    "四君子汤": ["Sijunzi Decoction", "Sijunzi Tang"],
    "四物汤": ["Siwu Decoction", "Siwu Tang"],
    "小柴胡汤": ["Xiaochaihu Decoction", "Xiaochaihu Tang", "Xiao Chai Hu"],
    "补中益气汤": ["Buzhong Yiqi Decoction", "Buzhong Yiqi Tang"],
    "桂枝汤": ["Guizhi Decoction", "Guizhi Tang"],
    "麻黄汤": ["Mahuang Decoction", "Mahuang Tang"],
    "八珍汤": ["Bazhen Decoction", "Bazhen Tang"],
    "黄连解毒汤": ["Huanglian Jiedu Decoction", "Huanglian Jiedu Tang"],
    "二陈汤": ["Erchen Decoction", "Erchen Tang"],
    "理中丸": ["Lizhong Wan", "Lizhong pill"],
    "六君子汤": ["Liujunzi Decoction", "Liujunzi Tang"],
    "当归补血汤": ["Danggui Buxue Decoction", "Danggui Buxue Tang",
                   "Danggui Buxue"],
    "寒喘祖帕颗粒": ["Hanchuan Zupa Granule", "Hanchuan Zupa"],
    "四味清肝汤": ["Siwei Qinggan Decoction", "Siwei Qinggan Tang"],
}

# 常见英文词若恰好是 HUGO 符号，直接剔除（避免摘要噪声）
_EN_WORDS = {"A", "I", "G", "T", "C", "AS", "AT", "IN", "OF", "ON", "OR",
             "TO", "UP", "BY", "FOR", "FROM", "THE", "AND", "ARE", "NOT",
             "HAS", "WAS", "WITH", "THAT", "THIS", "THAN", "OVER", "UNDER",
             "MANY", "SOME", "ALSO", "BETWEEN", "MOST", "LESS", "THROUGH",
             "WHICH", "WHERE", "WHEN", "AFTER", "BEFORE", "DURING", "AGAINST"}

_GENE_RE = re.compile(r"\b([A-Z][A-Z0-9]{1,24})\b")
_GENE_FULL = re.compile(r"[A-Z][A-Z0-9]{1,24}")
_YEAR_RE = re.compile(r"\b(19|20)\d{2}\b")

HGNC_SNAPSHOT = (Path(__file__).resolve().parents[3] / "data" / "reference" /
                 "hgnc_complete_set_2026-08-28.tsv")


def load_gene_dictionary(hgnc_path: Optional[Path] = None) -> set:
    """Load approved human protein-coding HGNC symbols.

    A dated official HGNC complete-set snapshot is preferred. Restricting the
    whitelist to approved ``gene with protein product`` entries matches this
    project's protein-target scope. The older project-local target dictionary
    remains an explicit fallback when the reference file is unavailable.
    """
    hgnc_path = Path(hgnc_path) if hgnc_path else HGNC_SNAPSHOT
    if hgnc_path.exists():
        genes = set()
        with hgnc_path.open(encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                symbol = (row.get("symbol") or "").strip()
                if (row.get("status") == "Approved" and
                        row.get("locus_type") == "gene with protein product" and
                        _GENE_FULL.fullmatch(symbol)):
                    genes.add(symbol)
        return genes - _EN_WORDS

    from .reproduce import PREF2GENE, PATHWAYS
    from .disease_focus import available_foci, resolve_focus
    genes = set()
    for gs in PATHWAYS.values():
        genes |= gs
    for theme in available_foci():
        try:
            genes |= set(resolve_focus(theme).get("targets", {}).keys())
        except Exception:  # noqa: BLE001  LLM 主题解析失败不影响词典
            pass
    genes |= set(PREF2GENE.values())
    genes -= _EN_WORDS
    return genes


def esearch_pmids(term: str, retmax: int = 5, sleep_s: float = 0.4) -> list:
    r = requests.get(f"{EUTILS}/esearch.fcgi",
                     params={"db": "pubmed", "term": term, "retmax": retmax,
                             "retmode": "json", "sort": "relevance"},
                     timeout=30)
    r.raise_for_status()
    time.sleep(sleep_s)
    return r.json().get("esearchresult", {}).get("idlist", [])


_MONTHS = {
    "jan": 1, "feb": 2, "mar": 3, "apr": 4, "may": 5, "jun": 6,
    "jul": 7, "aug": 8, "sep": 9, "oct": 10, "nov": 11, "dec": 12,
}


def _normalise_month(value: Optional[str]) -> Optional[int]:
    """Return an ISO-compatible month without guessing an unknown value."""
    if not value:
        return None
    value = value.strip()
    if value.isdigit() and 1 <= int(value) <= 12:
        return int(value)
    return _MONTHS.get(value[:3].lower())


def _element_date(node: Optional[ET.Element]) -> Optional[str]:
    """Extract the most precise available ISO-like date from a PubMed node."""
    if node is None:
        return None
    year = (node.findtext("Year") or "").strip()
    if not (year.isdigit() and len(year) == 4):
        medline = (node.findtext("MedlineDate") or "").strip()
        match = _YEAR_RE.search(medline)
        return match.group(0) if match else None
    month = _normalise_month(node.findtext("Month"))
    day = (node.findtext("Day") or "").strip()
    if month is None:
        return year
    if not (day.isdigit() and 1 <= int(day) <= 31):
        return f"{year}-{month:02d}"
    return f"{year}-{month:02d}-{int(day):02d}"


def _publication_date_fields(art: ET.Element) -> dict:
    """Apply one explicit policy for the temporal benchmark.

    The canonical date is the earliest PubMed ``ArticleDate`` (typically the
    article's online-first date), falling back to the journal issue date. PubMed
    search indexing may use a different date representation, so the benchmark
    retrieves one broad historical window and assigns records only after this
    common policy is applied.
    """
    article_dates = sorted(filter(None, (
        _element_date(node) for node in art.findall(".//Article/ArticleDate")
    )))
    issue_date = _element_date(art.find(".//Journal/JournalIssue/PubDate"))
    canonical = article_dates[0] if article_dates else issue_date
    return {
        "publication_date_policy": "first_article_date_else_journal_issue_date",
        "article_dates": article_dates,
        "journal_issue_date": issue_date,
        "canonical_publication_date": canonical,
        "canonical_publication_year": int(canonical[:4]) if canonical else None,
    }


def efetch_articles(pmids: List[str], sleep_s: float = 0.4) -> Dict[str, dict]:
    """Fetch PubMed XML with article-level provenance.

    A dated first-publication policy is used for the temporal split. Records
    without a parseable date remain in the raw audit trail but cannot enter a
    temporal window. DOI, journal and publication types are retained for
    auditable screening.
    """
    if not pmids:
        return {}
    r = requests.get(f"{EUTILS}/efetch.fcgi",
                     params={"db": "pubmed", "id": ",".join(pmids),
                             "rettype": "abstract", "retmode": "xml"},
                     timeout=30)
    r.raise_for_status()
    time.sleep(sleep_s)
    out: Dict[str, dict] = {}
    try:
        root = ET.fromstring(r.text)
    except ET.ParseError:
        return out
    for art in root.iter("PubmedArticle"):
        pmid = art.findtext(".//PMID")
        if not pmid:
            continue
        title = " ".join((art.findtext(".//ArticleTitle") or "").split())
        abstracts = ["".join(at.itertext())
                     for at in art.findall(".//Abstract/AbstractText")]
        doi = None
        for aid in art.findall(".//ArticleIdList/ArticleId"):
            if (aid.attrib.get("IdType") or "").lower() == "doi":
                doi = (aid.text or "").strip() or None
                break
        publication_types = sorted({
            " ".join("".join(pt.itertext()).split())
            for pt in art.findall(".//PublicationTypeList/PublicationType")
            if "".join(pt.itertext()).strip()
        })
        date_fields = _publication_date_fields(art)
        out[pmid] = {"title": title,
                     "abstract": " ".join(abstracts),
                     # ``year`` remains for backwards-compatible consumers;
                     # the canonical date fields are authoritative.
                     "year": date_fields["canonical_publication_year"],
                     "journal": " ".join(
                         (art.findtext(".//Journal/Title") or "").split()),
                     "doi": doi,
                     "publication_types": publication_types,
                     **date_fields}
    return out


def extract_genes(text: str, dictionary: set, min_count: int = 1) -> dict:
    """Extract dictionary genes verbatim from abstract text.

    Returns {gene: count}; genes with count < min_count are dropped.
    Case-sensitive; 'TNF-α' still matches 'TNF' because '-' is a word
    boundary.
    """
    counts: dict = {}
    for m in _GENE_RE.finditer(text):
        g = m.group(1)
        if g in dictionary:
            counts[g] = counts.get(g, 0) + 1
    return {g: c for g, c in sorted(counts.items()) if c >= min_count}


def extract_genes_llm(text: str, dictionary: set, llm=None,
                      min_count: int = 1) -> dict:
    """LLM-assisted gene extraction, whitelist-filtered.

    ``llm`` exposes ``extract_abstract_genes(text) -> {"genes": [{"gene",
    "assertion"}]}`` (or a dict with ``error``). A suggested gene is kept
    only if it appears verbatim in the abstract text AND matches HUGO
    symbol syntax — fabrication is structurally impossible. On LLM
    failure or an empty whitelist pass, falls back to the deterministic
    dictionary extraction.

    Returns {gene: {"count": n, "assertion": str}}.
    """
    def _dict_fallback() -> dict:
        return {g: {"count": c, "assertion": "dictionary"}
                for g, c in extract_genes(text, dictionary,
                                          min_count).items()}

    if llm is None:
        return _dict_fallback()
    try:
        resp = llm.extract_abstract_genes(text)
        items = resp.get("genes", []) if isinstance(resp, dict) else []
    except Exception:  # noqa: BLE001  LLM 失败退回字典法，绝不虚构
        return _dict_fallback()
    kept: dict = {}
    for item in items:
        gene = str(item.get("gene", "")).strip().upper()
        if not _GENE_FULL.fullmatch(gene):
            continue
        if gene in _EN_WORDS:
            continue
        if not re.search(rf"\b{re.escape(gene)}\b", text):
            continue  # 不在摘要逐字出现 → 丢弃
        count = len(re.findall(rf"\b{re.escape(gene)}\b", text))
        if count < min_count:
            continue
        kept[gene] = {"count": count,
                      "assertion": str(item.get("assertion", "mentioned"))}
    return kept or _dict_fallback()


def infer_pathways(genes: dict, min_hits: int = 2) -> list:
    from .reproduce import PATHWAYS
    out = []
    for pname, gs in PATHWAYS.items():
        hits = sorted(gs & set(genes))
        if len(hits) >= min_hits:
            out.append(pname)
    return out


def split_by_year(records: Dict[str, dict], cutoff: int = 2020
                  ) -> Tuple[Dict[str, dict], Dict[str, dict]]:
    """(train, test) by publication year relative to ``cutoff``.

    Articles without a parseable year are excluded from both windows.  Assigning
    an unknown date to development can make a later temporal comparison appear
    disjoint while silently changing the information set.  The v2 Paper A
    builder retains those records in its raw audit trail with an explicit reason.
    """
    train = {p: r for p, r in records.items()
             if (r.get("year") is not None) and (r.get("year") or 0) <= cutoff}
    test = {p: r for p, r in records.items()
            if (r.get("year") is not None) and (r.get("year") or 0) > cutoff}
    return train, test


def mine_formula(name: str, queries: List[str], retmax: int = 5,
                 min_count: int = 2, min_pathway_hits: int = 2,
                 cutoff: int = 2020, sleep_s: float = 0.4,
                 llm=None) -> Optional[dict]:
    """Mine one formula's PubMed benchmark with temporal split.

    Returns a benchmark dict (``lit_genes`` / ``lit_pathways`` /
    ``neg_pathways`` / ``benchmark_source`` incl. ``temporal_split``),
    or None when PubMed returns no PMIDs (honest skip, never guessed).
    """
    pmids: List[str] = []
    for q in queries:
        ids = esearch_pmids(f'"{q}" AND network pharmacology', retmax, sleep_s)
        pmids += ids
        if ids:
            break
    if not pmids:
        return None
    pmids = list(dict.fromkeys(pmids))[:retmax]
    records = efetch_articles(pmids, sleep_s)
    if not records:
        return None
    dictionary = load_gene_dictionary()
    total_text = " ".join(f"{r['title']} {r['abstract']}"
                          for r in records.values())
    if llm is None:
        genes = {g: {"count": c, "assertion": "dictionary"}
                 for g, c in extract_genes(total_text, dictionary,
                                           min_count).items()}
    else:
        genes = extract_genes_llm(total_text, dictionary, llm, min_count)
    lit_genes = sorted(genes)
    pathways = infer_pathways(genes, min_hits=min_pathway_hits)

    def _window_genes(recs: Dict[str, dict]) -> List[str]:
        text = " ".join(f"{r['title']} {r['abstract']}"
                        for r in recs.values())
        if llm is None:
            g = extract_genes(text, dictionary, min_count)
        else:
            g = extract_genes_llm(text, dictionary, llm, min_count)
        return sorted(g)

    train_recs, test_recs = split_by_year(records, cutoff)
    from .reproduce import PATHWAYS
    return {
        "lit_genes": lit_genes,
        "lit_pathways": pathways,
        "neg_pathways": sorted(set(PATHWAYS) - set(pathways)),
        "benchmark_source": {
            "query": queries[0],
            "pmids": pmids,
            "n_abstracts": len(records),
            "extracted_at": date.today().isoformat(),
            "method": ("PubMed E-utilities esearch+efetch；字典匹配（仅收录"
                       "摘要中逐字出现的基因）"
                       + ("；LLM 辅助抽取（白名单过滤，需逐字出现在摘要）"
                          if llm else "")),
            "temporal_split": {
                "cutoff": cutoff,
                "train_pmids": sorted(train_recs),
                "test_pmids": sorted(test_recs),
                "train_years": {p: r.get("year")
                                for p, r in sorted(train_recs.items())},
                "test_years": {p: r.get("year")
                               for p, r in sorted(test_recs.items())},
                "genes_train": _window_genes(train_recs),
                "genes_test": _window_genes(test_recs),
            },
        },
    }


def update_formula_library(entries: Dict[str, dict],
                           lib_path: Optional[Path] = None) -> int:
    """Write mined benchmarks into configs/formula_library.json."""
    if lib_path is None:
        lib_path = Path(__file__).resolve().parents[3] / "configs" / \
            "formula_library.json"
    lib = json.loads(lib_path.read_text(encoding="utf-8"))
    n = 0
    for name, bench in entries.items():
        if name not in lib["formulas"]:
            continue
        lib["formulas"][name].update(bench)
        n += 1
    lib_path.write_text(json.dumps(lib, ensure_ascii=False, indent=1),
                        encoding="utf-8")
    return n
