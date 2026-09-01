"""Three-tier evidence behavior-analysis engine.

Given one or more formulas (each: herbs -> documented compounds,
literature-anchored pathways and gene signatures), resolves
compound-target evidence under the ethno_evidence evidence architecture:

- measured tier:     ChEMBL measured activities, pchembl_value >= 5.0
                     (~10 µM), potency-ranked top 30, human protein targets
- weak_measured:     4.0 <= pchembl_value < 5.0 (~10-100 µM), sensitivity
- measured (TTD):    TTD local TSV as an independent measured source,
                     IC50/Ki/EC50 converted to µM, same bucketing function
- predicted tier:    SwissTargetPrediction structural-similarity prediction
                     (SEA principle); reported separately, never conflated

Wording is restricted to "consistent with / partially recapitulate";
no novelty or mechanism-confirmation claims are made.
"""
from __future__ import annotations

import json
import re
import threading
import time
from pathlib import Path
from typing import Dict, List, Optional, Set

import requests

from ..core.potency import PotencyRecord, bucket_by_potency, pchembl_to_um
from ..data.chembl_client import ChEMBLClient
from ..data.pubchem_client import PubChemClient
from ..data.swiss_client import SwissTargetPredictionClient
from ..data.ttd_client import TTDClient

# ---------------------------------------------------------------------------
# ChEMBL target preferred-name -> HUGO gene symbol normalization (curated,
# deterministic; does not add evidence, only re-labels).
# ---------------------------------------------------------------------------
PREF2GENE = {
    "Serine/threonine-protein kinase mTOR": "MTOR",
    "RAC-alpha serine/threonine-protein kinase": "AKT1",
    "NF-kappa-B p65 subunit": "RELA",
    "Transcription factor p65": "RELA",
    "NF-kappa-B p105 subunit": "NFKB1",
    "Nuclear factor NF-kappa-B p100 subunit": "NFKB2",
    "NF-kappa-B inhibitor alpha": "NFKBIA",
    "Inhibitor of nuclear factor kappa-B kinase subunit beta": "IKBKB",
    "Mitogen-activated protein kinase 14": "MAPK14",
    "Mitogen-activated protein kinase 1": "MAPK1",
    "Mitogen-activated protein kinase 3": "MAPK3",
    "Mitogen-activated protein kinase 8": "MAPK8",
    "Mitogen-activated protein kinase 9": "MAPK9",
    "Mitogen-activated protein kinase 10": "MAPK10",
    "Dual specificity mitogen-activated protein kinase kinase 1": "MAP2K1",
    "Dual specificity mitogen-activated protein kinase kinase 2": "MAP2K2",
    "Phosphatidylinositol 4,5-bisphosphate 3-kinase catalytic subunit alpha isoform": "PIK3CA",
    "Phosphatidylinositol 3-kinase regulatory subunit alpha": "PIK3R1",
    "3-phosphoinositide-dependent protein kinase 1": "PDPK1",
    "Ribosomal protein S6 kinase beta-1": "RPS6KB1",
    "Glycogen synthase kinase-3 beta": "GSK3B",
    "Eukaryotic translation initiation factor 4E-binding protein 1": "EIF4EBP1",
    "Serine/threonine-protein kinase ULK1": "ULK1",
    "Microtubule-associated proteins 1A/1B light chain 3B": "MAP1LC3B",
    "Sequestosome-1": "SQSTM1",
    "Beclin-1": "BECN1",
    "Autophagy-related protein 5": "ATG5",
    "Tumor necrosis factor": "TNF",
    "Interleukin-6": "IL6",
    "Interleukin-1 beta": "IL1B",
    "Interleukin-17A": "IL17A",
    "Interleukin-17F": "IL17F",
    "Interleukin-17 receptor A": "IL17RA",
    "Interleukin-17 receptor C": "IL17RC",
    "Interleukin-8": "CXCL8",
    "C-C motif chemokine 2": "CCL2",
    "Epidermal growth factor receptor": "EGFR",
    "Serine/threonine-protein kinase B-raf": "BRAF",
    "RAF proto-oncogene serine/threonine-protein kinase": "RAF1",
    "GTPase HRas": "HRAS",
    "GTPase KRas": "KRAS",
    "Prostaglandin G/H synthase 2": "PTGS2",
    "Transient receptor potential cation channel subfamily V member 1": "TRPV1",
    "5-hydroxytryptamine receptor 3A": "HTR3A",
    "Nuclear factor erythroid 2-related factor 2": "NFE2L2",
    "Cannabinoid receptor 2": "CNR2",
    "Arachidonate 5-lipoxygenase": "ALOX5",
    "Polyunsaturated fatty acid 5-lipoxygenase": "ALOX5",
    "Polyunsaturated fatty acid lipoxygenase ALOX15": "ALOX15",
    "Nitric oxide synthase, inducible": "NOS2",
    "Microtubule-associated protein 1 light chain 3 beta": "MAP1LC3B",
    "Microtubule-associated protein 1A/1B light chain 3B": "MAP1LC3B",
    "MAP kinase p38 alpha": "MAPK14",
    "Mitogen-activated protein kinase p38 alpha": "MAPK14",
    "RAC-beta serine/threonine-protein kinase": "AKT2",
    "Phosphatidylinositol 4,5-bisphosphate 3-kinase catalytic subunit beta isoform": "PIK3CB",
    "Dual specificity mitogen-activated protein kinase kinase 1/2": "MAP2K1",
    "Cellular tumor antigen p53": "TP53",
    "Cyclin-dependent kinase 1": "CDK1",
    "Cyclin-dependent kinase 1/cyclin B": "CDK1",
    "Cyclin-dependent kinase 2": "CDK2",
    "Cyclin-dependent kinase 5/CDK5 activator 1": "CDK5",
    "Cyclin-dependent kinase 6": "CDK6",
    "Casein kinase II subunit alpha": "CSNK2A1",
    "Casein kinase 2": "CSNK2A1",
    "Hypoxia-inducible factor 1-alpha": "HIF1A",
    "Aryl hydrocarbon receptor": "AHR",
    "Estrogen receptor": "ESR1",
    "Estrogen receptor beta": "ESR2",
    "Androgen receptor": "AR",
    "Thyrotropin receptor": "TSHR",
    "Tyrosine-protein kinase Fyn": "FYN",
    "Tyrosine-protein kinase SYK": "SYK",
    "Receptor-type tyrosine-protein kinase FLT3": "FLT3",
    "Receptor-type tyrosine-protein kinase C": "PTPRC",
    "Serine/threonine-protein kinase pim-1": "PIM1",
    "Serine/threonine-protein kinase NEK2": "NEK2",
    "Serine/threonine-protein kinase PLK1": "PLK1",
    "Serine/threonine-protein kinase Chk1": "CHEK1",
    "Serine/threonine-protein kinase Aurora-A": "AURKA",
    "Serine/threonine-protein kinase PAK 1": "PAK1",
    "Serine/threonine-protein kinase Sgk1": "SGK1",
    "Pyruvate kinase PKM": "PKM",
    "Ectonucleotide pyrophosphatase/phosphodiesterase family member 1": "ENPP1",
    "G-protein coupled receptor 35": "GPR35",
    "Inositol hexakisphosphate kinase 2": "IP6K2",
    "G-protein coupled bile acid receptor 1": "GPBAR1",
    "Nuclear receptor ROR-gamma": "RORC",
    "Tyrosine-protein phosphatase non-receptor type 1": "PTPN1",
    "Tyrosine-protein phosphatase non-receptor type 2": "PTPN2",
    "M-phase inducer phosphatase 2": "CDC25B",
    "Tissue factor": "F3",
    "Alpha-2A adrenergic receptor": "ADRA2A",
    "Alpha-2C adrenergic receptor": "ADRA2C",
    "Neuromedin-U receptor 2": "NMUR2",
    "GABA-A receptor; anion channel": "GABRA1",
    "Glycogen synthase kinase-3": "GSK3B",
    "Polyunsaturated fatty acid lipoxygenase ALOX12": "ALOX12",
    "Vascular endothelial growth factor receptor 2": "KDR",
    "Vascular endothelial growth factor receptor 1": "FLT1",
    "Fibroblast growth factor receptor 1": "FGFR1",
    "Tyrosine-protein kinase JAK2": "JAK2",
    "Tyrosine-protein kinase ABL1": "ABL1",
    "Bile acid receptor": "NR1H4",
    "Cannabinoid receptor 1": "CNR1",
    "Histone deacetylase 1": "HDAC1",
    "Prostaglandin G/H synthase 1": "PTGS1",
    "Thioredoxin reductase 1, cytoplasmic": "TXNRD1",
    "Transient receptor potential cation channel subfamily A member 1": "TRPA1",
    "ATP-binding cassette sub-family B member 1": "ABCB1",
    "Multidrug resistance protein 1": "ABCB1",
    "Acetylcholinesterase": "ACHE",
    "Butyrylcholinesterase": "BCHE",
    "5-hydroxytryptamine receptor 2A": "HTR2A",
    "D(2) dopamine receptor": "DRD2",
    "Sodium channel protein type 5 subunit alpha": "SCN5A",
    "Carbonic anhydrase 2": "CA2",
    "Carbonic anhydrase 1": "CA1",
    "Aldo-keto reductase family 1 member B1": "AKR1B1",
    "Aldo-keto reductase family 1 member A1": "AKR1A1",
    "Aldo-keto reductase family 7 member A3": "AKR7A3",
    "Aminopeptidase N": "ANPEP",
    "Neprilysin": "MME",
    "Sorbitol dehydrogenase": "SORD",
    "Glutathione reductase": "GSR",
    "Glutathione reductase, mitochondrial": "GSR",
    "cAMP-specific 3,5-cyclic phosphodiesterase 4D": "PDE4D",
    "Phosphodiesterase 5A": "PDE5A",
    "Adenosine A1 receptor": "ADORA1",
    "Adenosine A2a receptor": "ADORA2A",
    "Beta-2 adrenergic receptor": "ADRB2",
    "Histamine H1 receptor": "HRH1",
    "Kappa-type opioid receptor": "OPRK1",
    "Mu-type opioid receptor": "OPRM1",
    "Peroxisome proliferator-activated receptor alpha": "PPARA",
    "Liver X receptor alpha": "NR1H3",
    "Estrogen-related receptor alpha": "ESRRA",
    "Steroid 5-alpha-reductase 1": "SRD5A1",
    "Cytochrome P450 1A2": "CYP1A2",
    "Cytochrome P450 2C9": "CYP2C9",
    "Cytochrome P450 2D6": "CYP2D6",
    "Cytochrome P450 3A4": "CYP3A4",
    "Cytochrome P450 1A1": "CYP1A1",
    "Cytochrome P450 2E1": "CYP2E1",
    "Cytochrome P450 2C19": "CYP2C19",
    "UDP-glucuronosyltransferase 1A1": "UGT1A1",
    "UDP-glucuronosyltransferase 1A4": "UGT1A4",
    "Glutathione S-transferase P": "GSTP1",
    "Alpha-1-acid glycoprotein 1": "ORM1",
    "Serum albumin": "ALB",
    "Transthyretin": "TTR",
    "Heat shock protein HSP 90-alpha": "HSP90AA1",
    "Heat shock protein HSP 90-beta": "HSP90AB1",
    "Peptidyl-prolyl cis-trans isomerase A": "PPIA",
    "Beta-secretase 1": "BACE1",
    "Fatty acid synthase": "FASN",
    "5'-AMP-activated protein kinase catalytic subunit alpha-1": "PRKAA1",
    "5'-AMP-activated protein kinase catalytic subunit alpha-2": "PRKAA2",
    "HMG-CoA reductase": "HMGCR",
    "Dipeptidyl peptidase 4": "DPP4",
    "Angiotensin-converting enzyme": "ACE",
    "Renin": "REN",
    "Diacylglycerol O-acyltransferase 1": "DGAT1",
    "Fatty-acid amide hydrolase 1": "FAAH",
    "Aldehyde dehydrogenase 1A1": "ALDH1A1",
    "Thyroid hormone receptor beta": "THRB",
    "Olfactory receptor 5K1": "OR5K1",
    "3-oxo-5-alpha-steroid 4-dehydrogenase 1": "SRD5A1",
    "11-beta-hydroxysteroid dehydrogenase 1": "HSD11B1",
    "11-beta-hydroxysteroid dehydrogenase type 2": "HSD11B2",
    "15-hydroxyprostaglandin dehydrogenase [NAD(+)]": "HPGD",
    "17-beta-hydroxysteroid dehydrogenase type 1": "HSD17B1",
    "17-beta-hydroxysteroid dehydrogenase type 2": "HSD17B2",
    "3-hydroxyacyl-CoA dehydrogenase type-2": "HADH",
    "72 kDa type IV collagenase": "MMP2",
    "ATP-dependent DNA helicase Q1": "RECQL",
    "Aldo-keto reductase family 1 member B10": "AKR1B10",
    "Alpha-(1,3)-fucosyltransferase 7": "FUT7",
    "Alpha-synuclein": "SNCA",
    "Amine oxidase [flavin-containing] A": "MAOA",
    "Aromatase": "CYP19A1",
    "C-C chemokine receptor type 4": "CCR4",
    "CDGSH iron-sulfur domain-containing protein 1": "CISD1",
    "Carbonic anhydrase 12": "CA12",
    "Carbonic anhydrase 14": "CA14",
    "Carbonic anhydrase 3": "CA3",
    "Carbonic anhydrase 4": "CA4",
    "Carbonic anhydrase 5A, mitochondrial": "CA5A",
    "Carbonic anhydrase 5B, mitochondrial": "CA5B",
    "Carbonic anhydrase 6": "CA6",
    "Carbonic anhydrase 7": "CA7",
    "Carbonic anhydrase 9": "CA9",
    "Chymotrypsinogen B": "CTRB1",
    "Chymotrypsin-C": "CTRC",
    "Cytochrome P450 1B1": "CYP1B1",
    "DNA polymerase beta": "POLB",
    "DNA repair nuclease/redox regulator APEX1": "APEX1",
    "High mobility group protein B1": "HMGB1",
    "Induced myeloid leukemia cell differentiation protein Mcl-1": "MCL1",
    "Interstitial collagenase": "MMP1",
    "Lysine-specific demethylase 4E": "KDM4E",
    "Lysine-specific histone demethylase 1A": "KDM1A",
    "Lysosomal alpha-glucosidase": "GAA",
    "Macrophage metalloelastase": "MMP12",
    "Menin/Histone-lysine N-methyltransferase MLL": "MEN1",
    "Microtubule-associated protein tau": "MAPT",
    "Multidrug resistance-associated protein 1": "ABCC1",
    "Myeloperoxidase": "MPO",
    "NADPH oxidase 4": "NOX4",
    "Poly [ADP-ribose] polymerase 1": "PARP1",
    "Poly [ADP-ribose] polymerase tankyrase-2": "TNKS2",
    "Polypeptide N-acetylgalactosaminyltransferase 2": "GALNT2",
    "Prelamin-A/C": "LMNA",
    "Protein deacetylase HDAC6": "HDAC6",
    "Prothrombin": "F2",
    "Tyrosinase": "TYR",
    "Tyrosyl-DNA phosphodiesterase 1": "TDP1",
    "Xanthine dehydrogenase/oxidase": "XDH",
    "Hypoxia-inducible factor 1-alpha": "HIF1A",
    "Integrin alpha-5": "ITGA5",
    "Integrin beta-1": "ITGB1",
    "Integrin alpha-V": "ITGAV",
    "Vascular endothelial growth factor A": "VEGFA",
    "Focal adhesion kinase 1": "PTK2",
    "Proto-oncogene tyrosine-protein kinase Src": "SRC",
    "Fibronectin": "FN1",
    "Talin-1": "TLN1",
    "NACHT, LRR and PYD domains-containing protein 3": "NLRP3",
    "Caspase-1": "CASP1",
    "Interleukin-18": "IL18",
    "Gasdermin-D": "GSDMD",
    "Apoptosis-associated speck-like protein containing a CARD": "PYCARD",
    "Toll-like receptor 4": "TLR4",
    "L-lactate dehydrogenase A chain": "LDHA",
    "Erythropoietin": "EPO",
    "Metalloproteinase inhibitor 1": "TIMP1",
    "Solute carrier family 2, facilitated glucose transporter member 1": "SLC2A1",
    "Solute carrier family 2, facilitated glucose transporter member 4": "SLC2A4",
    "Pyruvate dehydrogenase (acetyl-transferring) kinase isozyme 1, mitochondrial": "PDK1",
}

# ChEMBL 'target' labels that are cell lines or non-protein assay labels;
# excluded after the organism filter to keep the protein-target view honest.
NOISE_TARGETS = {
    "A549", "ANN-1", "BGC-823", "BXPC-3", "Calu-6", "DU-145", "HFF", "HL-60",
    "HT-1080", "HeLa", "Huh-7", "K562", "KM-20L2", "MCF7", "MOLM-13",
    "MOLM-14", "NCI-H460", "NIH3T3", "P388", "RAW", "SF-268", "THP-1",
    "Vero", "WiDr",
    "Unchecked", "No relevant target", "Cell membrane", "ADMET", "Unknown",
    "Monocyte", "microRNA 21",
}

_IL_RE = re.compile(r"^Interleukin-(\d+)([A-Z]?)$")


def normalize(pref_name: str) -> Optional[str]:
    """Map a ChEMBL target preferred name to a HUGO gene symbol; None if unmapped."""
    if not pref_name:
        return None
    if pref_name in PREF2GENE:
        return PREF2GENE[pref_name]
    m = _IL_RE.match(pref_name)
    if m:
        return "IL" + m.group(1) + (m.group(2) or "")
    return None


# ---------------------------------------------------------------------------
# Pathway gene sets (KEGG-informed core members; used only for the
# pathway-consistency check, not claimed as exhaustive).
# ---------------------------------------------------------------------------
PATHWAYS = {
    "NF-kB signaling": {"RELA", "NFKB1", "NFKB2", "NFKBIA", "IKBKB", "CHUK",
                        "TNF", "IL1B", "IL6", "CCL2", "CXCL8"},
    "MAPK signaling": {"MAPK1", "MAPK3", "MAPK8", "MAPK9", "MAPK10", "MAPK14",
                       "MAP2K1", "MAP2K2", "BRAF", "RAF1", "HRAS", "KRAS",
                       "EGFR", "TNF"},
    "PI3K-Akt-mTOR": {"PIK3CA", "PIK3R1", "AKT1", "MTOR", "PDPK1", "RPS6KB1",
                      "EIF4EBP1", "GSK3B", "PTEN", "TSC2", "HIF1A", "TP53",
                      "SGK1"},
    "Autophagy": {"ULK1", "MAP1LC3B", "SQSTM1", "BECN1", "ATG5", "ATG7",
                  "ATG12", "LAMP1"},
    "IL-17 signaling": {"IL17A", "IL17F", "IL17RA", "IL17RC", "IL6", "CXCL8",
                        "CCL2", "TNF", "MAPK1", "MAPK3", "MAPK8", "NFKB1",
                        "RELA", "IKBKB"},
    "Cytokine (predicted)": {"TNF", "IL1B", "IL6", "CXCL8"},
    "HIF-1 signaling": {"HIF1A", "VEGFA", "FLT1", "KDR", "PDK1", "LDHA",
                        "SLC2A1", "NOS2", "EPO", "TIMP1"},
    "Integrin signaling": {"ITGA5", "ITGB1", "ITGAV", "PTK2", "SRC", "FN1",
                           "TLN1", "VEGFA"},
    "NLRP3 inflammasome": {"NLRP3", "CASP1", "PYCARD", "IL1B", "IL18",
                           "GSDMD", "TLR4", "NFKB1", "RELA"},
}


def jaccard(a: Set[str], b: Set[str]) -> float:
    if not a or not b:
        return 0.0
    return round(len(a & b) / len(a | b), 4)


def compute_metrics(meta: dict, derived: Set[str]) -> dict:
    """Recapture + Jaccard + pathway-consistency for one gene set.

    Also reports negative-control pathway hits (``neg_pathways`` in the
    formula metadata): pathways the formula's literature did NOT validate,
    used as a specificity reference. Negative hits are reported as-is; a
    nonzero value does not disprove the mechanism, it only lowers the
    specificity signal of that tier.
    """
    lit = set(meta["lit_genes"])
    recap = sorted(lit & derived)
    ph = {}
    for pname, genes in PATHWAYS.items():
        if pname in meta["lit_pathways"]:
            hits = sorted(derived & genes)
            ph[pname] = {"hit_count": len(hits), "hits": hits}
    neg_ph = {}
    for pname, genes in PATHWAYS.items():
        if pname in meta.get("neg_pathways", []):
            hits = sorted(derived & genes)
            neg_ph[pname] = {"hit_count": len(hits), "hits": hits}
    return {
        "derived_gene_count": len(derived),
        "recaptured_genes": recap,
        # 无文献基准（lit_genes 为空）时复现率=None，前端显示 N/A，
        # 避免 0/0 误读为"复现失败"
        "recapture_rate": (round(len(recap) / len(lit), 4) if lit else None),
        "jaccard": jaccard(lit, derived),
        "pathway_hits": ph,
        "negative_pathway_hits": neg_ph,
        "negative_total_hits": sum(v["hit_count"] for v in neg_ph.values()),
    }


def resolve_with_scores(client, name: str, herb: Optional[str] = None,
                        organism: str = "Homo sapiens",
                        fetch_limit: int = 300) -> dict:
    """One-fetch, locally bucketed ChEMBL target resolution.

    Strong tier: potency <= 10 uM (pchembl_value >= 5.0), top 30 by potency.
    Weak tier: 10-100 uM (pchembl_value 4.0-5.0). Bucketing shares the
    unified criterion in core.potency (identical to the TTD path); the two
    tiers are reported separately and never merged.

    When ``herb`` is given, the offline snapshot (chembl_measured.json) is
    checked first: hits resolve instantly with the same bucketing already
    applied. Live ChEMBL is only the fallback for snapshot misses.
    """
    if herb is not None:
        from ..data.chembl_snapshot import resolve_from_snapshot
        hit = resolve_from_snapshot(herb, name)
        if hit:
            return hit
    try:
        mols = client.search_molecules(name, limit=1)
        if not mols:
            return {"chembl_id": None, "strong": [], "weak": []}
        chembl_id = mols[0]["molecule_chembl_id"]
        acts = client.target_activities(chembl_id, limit=fetch_limit,
                                        organism=organism,
                                        order_by="-pchembl_value")
        records = [
            PotencyRecord(
                target=a["target_pref_name"],
                potency_um=pchembl_to_um(a.get("pchembl_value") or 0.0),
                activity_type=a.get("standard_type") or "",
                raw=str(a.get("pchembl_value")),
            )
            for a in acts
            if a.get("target_pref_name") and a.get("pchembl_value") is not None
        ]
        strong, weak = bucket_by_potency(records)
        time.sleep(0.35)
        return {"chembl_id": chembl_id, "strong": strong, "weak": weak}
    except requests.RequestException:
        return {"chembl_id": None, "strong": [], "weak": [], "error": "api_failure"}


def _filter_noise(entries: List[dict]) -> tuple:
    """Remove noise labels from compound entries; return counts + entries."""
    n_noise = 0
    n_unmapped = 0
    for e in entries:
        for tier_key in ("strong", "weak"):
            kept = []
            for t in e.get(tier_key, []):
                if t in NOISE_TARGETS:
                    n_noise += 1
                else:
                    kept.append(t)
            e[tier_key] = kept
    return n_noise, n_unmapped


def run_reproduction(formulas: Dict[str, dict], with_stp: bool = False,
                     with_herb: bool = False, with_etcm: bool = False,
                     out_path: Optional[Path] = None,
                     verbose: bool = True,
                     stp_workers: int = 4) -> dict:
    """Run three-tier evidence behavior analysis over formulas.

    ``formulas``: dict of {formula_id: {"zh_name", "pmid", "doi",
    "indication", "herbs": {herb: [compounds]}, "lit_pathways": [...],
    "lit_genes": [...]}} — same schema as evidence.benchmarks.FORMULAS.

    Returns the full structured result dict; writes it to ``out_path`` if
    given. The predicted tier (STP) runs only when ``with_stp`` is True.
    When enabled, STP submissions run in a background thread in parallel
    with the measured/TTD computation (pipeline overlap), then results
    are joined — so the predicted tier adds almost no serial wall time
    on top of the measured tier. ``stp_workers`` controls the number of
    concurrent STP submissions (default 4; be mindful of server rate
    limits).
    """
    client = ChEMBLClient()
    out = {
        "note": "Three-tier evidence behavior analysis under the ethno_evidence "
                "evidence architecture. Filters: target_organism=Homo sapiens, "
                "pchembl_value >= 5.0 (~10 uM), potency-ranked top-30, "
                "cell-line/non-protein labels excluded; TTD (v10.1.01 local TSV) "
                "as independent measured source; predicted tier = SwissTarget"
                "Prediction structural-similarity prediction, reported separately. "
                "Wording limited to consistent-with / partially-recapitulate. "
                "Literature-anchored formula examples are behavioral illustrations, "
                "not independent ground truth or predictive-performance benchmarks.",
        "filters": {"organism": "Homo sapiens"},
        "formulas": {},
    }
    for key, meta in formulas.items():
        if verbose:
            print(f"\n=== {key} ({meta.get('zh_name', '')}) PMID {meta.get('pmid', '')} ===")
        all_compounds = [c for comps in meta["herbs"].values() for c in comps]
        # --- 预测层预提交：SMILES 并发解析 + STP 后台线程，与实测层并行 ---
        stp_thread = None
        stp_box: dict = {}
        if with_stp:
            from concurrent.futures import ThreadPoolExecutor
            pubchem = PubChemClient()

            def _smi(cname: str):
                return cname, pubchem.get_smiles(cname)

            name_smiles = {}
            with ThreadPoolExecutor(max_workers=6) as ex:
                for cname, smi in ex.map(_smi, set(all_compounds)):
                    if smi:
                        name_smiles[cname] = smi

            def _stp_run():
                stp_box["res"] = SwissTargetPredictionClient().predict_many(
                    name_smiles, workers=stp_workers)

            stp_thread = threading.Thread(target=_stp_run, daemon=True)
            stp_thread.start()
            if verbose:
                print(f"[STP      ] 后台预提交 {len(name_smiles)} 个化合物 "
                      f"（{stp_workers} 并发），与实测层并行…")
        # --- compound -> targets, bucketed by potency (strong >=5.0, weak 4.0-5.0) ---
        compound_entries = []
        for herb, compounds in meta["herbs"].items():
            for cname in compounds:
                entry = {"herb": herb, "compound": cname,
                         **resolve_with_scores(client, cname, herb=herb)}
                compound_entries.append(entry)
        # --- aggregate genes per evidence tier ---
        derived_strong: Set[str] = set()
        derived_weak_only: Set[str] = set()
        n_resolved = sum(1 for e in compound_entries if e.get("chembl_id"))
        n_failed = sum(1 for e in compound_entries if e.get("error") == "api_failure")
        n_unmapped = 0
        n_noise, _ = _filter_noise(compound_entries)
        for e in compound_entries:
            for tier_key in ("strong", "weak"):
                for pref in e.get(tier_key, []):
                    g = normalize(pref)
                    if g:
                        if tier_key == "strong":
                            derived_strong.add(g)
                        else:
                            derived_weak_only.add(g)
                    else:
                        n_unmapped += 1
        derived_all = derived_strong | derived_weak_only
        lit = set(meta["lit_genes"])
        m_strong = compute_metrics(meta, derived_strong)
        m_all = compute_metrics(meta, derived_all)
        # --- TTD path: same pipeline, local TSV, unified potency bucketing ---
        # 本地 TTD TSV（v10.1.01）缺失时如实降级：TTD 层标记不可用并记录
        # 具体缺失文件，不虚构任何解析结果，其余实测层继续运行——与
        # HERB/ETCM 层"本地未下载 TSV 时如实返回 0 解析"的语义保持一致。
        ttd_res = {}
        ttd_strong: Set[str] = set()
        ttd_weak_only: Set[str] = set()
        ttd_resolved = 0
        ttd_error = None
        try:
            ttd_res = TTDClient().resolve_compound_targets(all_compounds)
        except FileNotFoundError as exc:
            ttd_error = f"ttd_data_unavailable: {exc}"
            if verbose:
                print(f"[TTD      ] 本地 TTD TSV 缺失，TTD 层不可用 "
                      f"（如实记录，不虚构结果）：{exc}")
        else:
            for cname, e in ttd_res.items():
                if e.get("ttd_id"):
                    ttd_resolved += 1
                ttd_strong.update(e.get("strong", []))
                ttd_weak_only.update(e.get("weak", []))
        ttd_all = ttd_strong | ttd_weak_only
        mt_s = compute_metrics(meta, ttd_strong)
        mt_a = compute_metrics(meta, ttd_all)
        # --- HERB path: reference-mining tier (文献整合数据, 与实测严格分档) ---
        herb = {"enabled": with_herb,
                "evidence_level": "reference_mining",
                "compounds_queried": len(all_compounds)}
        if with_herb:
            from ..data.herb_client import HERBClient
            herb_res = HERBClient().resolve_compound_targets(all_compounds)
            herb_genes: Set[str] = set()
            herb_resolved = 0
            for cname, e in herb_res.items():
                if e.get("targets"):
                    herb_resolved += 1
                    herb_genes.update(
                        g for g in e["targets"]
                        if isinstance(g, str) and g and g not in NOISE_TARGETS)
            herb_lit = set(meta["lit_genes"])
            herb_recap = sorted(herb_lit & herb_genes)
            herb.update({
                "compounds_resolved_herb": herb_resolved,
                "gene_count": len(herb_genes),
                "genes": sorted(herb_genes),
                "recaptured_genes": herb_recap,
                "recapture_rate": (round(len(herb_recap) / len(herb_lit), 4)
                                   if herb_lit else None),
                "note": ("HERB 参考挖掘层（herb.ac.cn，文献整合数据）："
                         "与实测层严格分档，不参与复现率汇总；本地未下载 "
                         "HERB 导出 TSV 时如实返回 0 解析。"),
            })
            if verbose:
                print(f"[HERB      ] reference-mining resolved={herb_resolved}"
                      f"/{len(all_compounds)}  genes={len(herb_genes)}"
                      f"  recapture={len(herb_recap)}/{len(herb_lit)}")
        # --- ETCM path: herb-level reference mining (药材级预测靶点, 与实测严格分档) ---
        etcm = {"enabled": with_etcm, "evidence_level": "reference_mining",
                "herbs_queried": len(meta["herbs"])}
        if with_etcm:
            from ..data.etcm_client import ETCMClient
            etcm_res = ETCMClient().resolve_herb_targets(list(meta["herbs"]))
            etcm_genes: Set[str] = set()
            etcm_resolved = 0
            for hname, e in etcm_res.items():
                if e.get("targets"):
                    etcm_resolved += 1
                    etcm_genes.update(
                        g for g in e["targets"]
                        if isinstance(g, str) and g and g not in NOISE_TARGETS)
            etcm_lit = set(meta["lit_genes"])
            etcm_recap = sorted(etcm_lit & etcm_genes)
            etcm.update({
                "herbs_resolved_etcm": etcm_resolved,
                "gene_count": len(etcm_genes),
                "genes": sorted(etcm_genes),
                "recaptured_genes": etcm_recap,
                "recapture_rate": (round(len(etcm_recap) / len(etcm_lit), 4)
                                   if etcm_lit else None),
                "note": ("ETCM 药材级参考挖掘层（www.tcmip.cn/ETCM，计算预测候选靶点 "
                         "confidence>=0.80）：与实测层严格分档，不参与复现率汇总；"
                         "本地未下载 ETCM TSV 时如实返回 0 解析。"),
            })
            if verbose:
                print(f"[ETCM      ] herb-level resolved={etcm_resolved}"
                      f"/{len(meta['herbs'])}  genes={len(etcm_genes)}"
                      f"  recapture={len(etcm_recap)}/{len(etcm_lit)}")
        # --- STP path: structural-similarity prediction (predicted tier) ---
        stp = {"compounds_queried": len(all_compounds), "enabled": with_stp}
        if with_stp:
            stp_thread.join()  # 等待后台预提交完成（与实测层并行）
            stp_res = stp_box.get("res", {})
            stp_genes: Set[str] = set()
            stp_resolved = 0
            for cname, e in stp_res.items():
                if e.get("targets"):
                    stp_resolved += 1
                for t in e.get("targets", []):
                    g = t.get("common_name")
                    if g and g not in NOISE_TARGETS:
                        stp_genes.add(g)
            mstp = compute_metrics(meta, stp_genes)
            # 逐化合物预测靶点（chembl_id + 基因符号 + 概率）—— 校准
            # （uncertainty 模块，按 ChEMBL ID join 实测快照）与前瞻验证
            # （按基因符号 join 时间线切割后的新文献）共用，口径一致。
            # 概率为 None 时按 0.0 计入（校准曲线对这些样本给出最低置信）。
            compound_targets = {
                cname: {"targets": [
                    {"chembl_id": t.get("chembl_id"),
                     "gene": t.get("common_name"),
                     "probability": t.get("probability") or 0.0}
                    for t in e.get("targets", []) if t.get("chembl_id")]}
                for cname, e in stp_res.items()
                if any(t.get("chembl_id") for t in e.get("targets", []))
            }
            # Top30 文献重叠率：预测层高置信头部与预设文献集合的交集/30；
            # 文献集合并非完整金标准，因此不是准确率。该比例替代全量 recall，
            # 以避免大预测集撞小文献集导致覆盖率虚高。
            from .uncertainty import topk_predicted_genes, untested_fraction
            top30_genes = topk_predicted_genes(compound_targets, k=30)
            mstp["top30_precision"] = (
                round(len(set(top30_genes) & lit) / 30, 4) if lit else None)
            # 未测试率：预测靶点中从未在 ChEMBL 快照出现的比例（覆盖缺口，
            # 未命中≠阴性）
            from ..data.chembl_snapshot import all_measured_target_ids
            pred_ids = {t.get("chembl_id")
                        for e in stp_res.values()
                        for t in e.get("targets", [])
                        if t.get("chembl_id")}
            stp["untested_fraction"] = untested_fraction(
                pred_ids, all_measured_target_ids())
            stp.update({
                "smiles_obtained": len(name_smiles),
                "compounds_resolved_stp": stp_resolved,
                "predicted_gene_count": len(stp_genes),
                "strong_tier_only": mstp,
                "strong_plus_weak": mstp,
                "compound_targets": compound_targets,
            })
            if verbose:
                print(f"[STP      ] resolved compounds={stp_resolved}/{len(name_smiles)}  "
                      f"predicted genes={len(stp_genes)}")
                if lit:
                    print(f"[STP pred ] Top30文献重叠率 "
                          f"{len(set(top30_genes) & lit)}/30 "
                          f"({mstp['top30_precision']:.1%})  "
                          f"复现率(全量) {mstp['recapture_rate']:.1%}  "
                          f"Jaccard={mstp['jaccard']}")
                else:
                    print(f"[STP pred ] Top30文献重叠率 N/A（无文献集合）  "
                          f"预测基因 {len(stp_genes)}")
                print(f"[STP pred ] 未测试率 {stp['untested_fraction']:.1%}"
                      f"（预测靶点中 ChEMBL 未测过的比例；未命中≠阴性）")
        # --- report ---
        if verbose:
            def _rate(recap: list, lit: set) -> str:
                return (f"{len(recap)}/{len(lit)} "
                        f"({len(recap) / len(lit):.1%})" if lit else
                        f"{len(recap)}/0 (N/A 无文献基准)")
            print(f"herbs={len(meta['herbs'])}  compounds={len(compound_entries)}  "
                  f"resolved={n_resolved}  api_failures={n_failed}")
            print(f"strong genes (>=10 uM) = {len(derived_strong)}  "
                  f"(+weak 10-100 uM = {len(derived_all)})  "
                  f"unmapped = {n_unmapped}  noise removed = {n_noise}")
            print(f"[strong only ] recapture {_rate(m_strong['recaptured_genes'], lit)}"
                  f"  Jaccard={m_strong['jaccard']}")
            print(f"[strong+weak ] recapture {_rate(m_all['recaptured_genes'], lit)}"
                  f"  Jaccard={m_all['jaccard']}")
            if not ttd_error:
                print(f"[TTD       ] resolved compounds={ttd_resolved}/{len(all_compounds)}  "
                      f"strong genes={len(ttd_strong)}  (+weak={len(ttd_all)})")
                print(f"[TTD strong] recapture {_rate(mt_s['recaptured_genes'], lit)}"
                      f"  Jaccard={mt_s['jaccard']}")
        if meta.get("neg_pathways"):
            neg_keys = sorted(m_all["negative_pathway_hits"])
            if verbose:
                print("[negative-control] 未验证通路命中（特异性参考）:")
                for pname in neg_keys:
                    v = m_all["negative_pathway_hits"][pname]
                    print(f"  - {pname}: {v['hit_count']} "
                          f"({', '.join(v['hits']) if v['hits'] else 'none'})")
        out["formulas"][key] = {
            "zh_name": meta["zh_name"],
            "pmid": meta.get("pmid"),
            "doi": meta.get("doi"),
            "compounds_queried": len(compound_entries),
            "compounds_resolved_chembl": n_resolved,
            "api_failures": n_failed,
            "unmapped_target_labels": n_unmapped,
            "literature_gene_count": len(lit),
            "lit_genes": sorted(lit),
            "lit_pathways": meta.get("lit_pathways", []),
            "benchmark_source": meta.get("benchmark_source", {}) or {},
            "neg_pathways": meta.get("neg_pathways", []),
            "strong_tier_only": m_strong,
            "strong_plus_weak": m_all,
            "specificity": {
                "negative_control_note": (
                    "Pathways the formula's literature did not validate; "
                    "hits are reported as-is for specificity reference and do "
                    "not disprove the mechanism."
                ),
                "positive_pathway_total_hits": sum(
                    v["hit_count"] for v in m_all["pathway_hits"].values()
                ),
                "negative_pathway_total_hits": m_all["negative_total_hits"],
                "negative_pathway_detail": m_all["negative_pathway_hits"],
            },
            "compound_detail": compound_entries,
            "ttd": {
                "compounds_queried": len(all_compounds),
                "compounds_resolved_ttd": ttd_resolved,
                "literature_gene_count": len(lit),
                "error": ttd_error,
                "strong_tier_only": mt_s,
                "strong_plus_weak": mt_a,
                "compound_detail": {
                    cn: {"ttd_id": e.get("ttd_id"),
                         "n_records": len(e.get("records", [])),
                         "strong": e.get("strong", []),
                         "weak": e.get("weak", []),
                         "error": e.get("error")}
                    for cn, e in ttd_res.items()
                },
            },
            "herb": herb,
            "etcm": etcm,
            "stp": stp,
        }
    if out_path:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(json.dumps(out, ensure_ascii=False, indent=2),
                            encoding="utf-8")
        if verbose:
            print(f"\n[ok] results written to {out_path}")
    return out
