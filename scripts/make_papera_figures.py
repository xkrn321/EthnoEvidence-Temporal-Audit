"""Create publication figures for the temporally audited methodology manuscript.

The figures are deliberately descriptive: all values are read from the locked
v2 benchmark files and show evidence coverage or screening workload, not
predictive performance or biological efficacy.
"""
from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

ROOT = Path(__file__).resolve().parent.parent
BENCH = ROOT / "benchmarks" / "papera_temporal_v2"
# In the public repository, figures and machine-readable audits live beside
# the benchmark and are regenerated without controlled fieldwork files.
PACKAGE = ROOT
OUT = PACKAGE / "figures"
AUDITS = PACKAGE / "audits"
OUT.mkdir(parents=True, exist_ok=True)

BLUE = "#2F5D7E"
TEAL = "#2A9D8F"
ORANGE = "#E9A03B"
RED = "#B64A4A"
GREY = "#6B7280"
LIGHT = "#E9EEF2"
DARK = "#263238"

plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 8,
    "axes.titlesize": 10,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
})


def _save(fig, stem: str):
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.svg", bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.png", dpi=600, bbox_inches="tight")
    plt.close(fig)


def figure1_framework():
    # Keep this panel deliberately roomy: it is embedded at ~6 in width in
    # the manuscript, so long labels must be wrapped before rasterisation.
    fig, ax = plt.subplots(figsize=(7.1, 4.0))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    def box(x, y, w, h, title, lines, color=BLUE, dashed=False):
        patch = FancyBboxPatch(
            (x, y), w, h, boxstyle="round,pad=0.012,rounding_size=0.018",
            linewidth=1.2, edgecolor=color, facecolor="white",
            linestyle="--" if dashed else "-",
        )
        ax.add_patch(patch)
        # Every text artist is clipped to its own rounded patch.  This is the
        # final guard against glyphs escaping the border after DOCX/PDF
        # rescaling, while explicit line breaks preserve readability.
        title_lines = title.split("\n")
        title_artist = ax.text(x + w / 2, y + h - 0.035,
                               "\n".join(title_lines), ha="center", va="top",
                               fontsize=7.8 if len(title_lines) == 1 else 7.2,
                               fontweight="bold", color=color, linespacing=1.0,
                               clip_on=True)
        title_artist.set_clip_path(patch)
        # Body is kept inside a generous inset and uses short, pre-wrapped
        # lines rather than relying on renderer-dependent word wrapping.
        body_offset = 0.102 + 0.060 * (len(title_lines) - 1)
        body_artist = ax.text(x + 0.020, y + h - body_offset,
                              "\n".join(lines), ha="left", va="top",
                              fontsize=5.9, color=DARK, linespacing=1.05,
                              clip_on=True)
        body_artist.set_clip_path(patch)
        return patch

    def arrow(a, b, color=GREY, dashed=False):
        ax.add_patch(FancyArrowPatch(
            a, b, arrowstyle="-|>", mutation_scale=10, linewidth=1.1,
            color=color, linestyle="--" if dashed else "-"))

    box(0.03, 0.58, 0.18, 0.27, "Study input",
        ["Multicomponent formula", "Marker compounds", "Source provenance"], ORANGE)
    box(0.27, 0.58, 0.20, 0.27, "Calibration\nevidence",
        ["14-formula cohort", "PubMed / ChEMBL <= 2022", "STP comparator"], BLUE)
    box(0.53, 0.58, 0.18, 0.27, "Candidate\nrecord",
        ["Provenance", "Assay / potency", "Missing-data flags"], TEAL)
    box(0.77, 0.58, 0.20, 0.27, "Evidence-gated\nrank",
        ["Tier before score", "Deterministic tie-\nbreak", "Claim ceiling"], RED)
    arrow((0.21, 0.715), (0.27, 0.715))
    arrow((0.47, 0.715), (0.53, 0.715))
    arrow((0.71, 0.715), (0.77, 0.715))

    ax.plot([0.5, 0.5], [0.42, 0.47], color=RED, linewidth=1.5,
            linestyle="--")
    boundary_patch = FancyBboxPatch(
        (0.34, 0.47), 0.32, 0.07,
        boxstyle="round,pad=0.008,rounding_size=0.014",
        linewidth=1.0, edgecolor=RED, facecolor="white", linestyle="--")
    ax.add_patch(boundary_patch)
    boundary_text = ax.text(0.50, 0.505, "Frozen temporal boundary",
                            ha="center", va="center", color=RED,
                            fontsize=7.4, fontweight="bold", clip_on=True)
    boundary_text.set_clip_path(boundary_patch)
    box(0.03, 0.22, 0.40, 0.22, "Development / calibration",
        ["No post-cutoff labels used", "Locked cohort + baselines",
         "n = 593 formula-article records"],
        BLUE)
    box(0.57, 0.22, 0.40, 0.22, "Held-out literature audit",
        ["PubMed 2023-2025", "Blinded review + adjudication",
         "n = 339; future evidence"], TEAL, dashed=True)
    arrow((0.43, 0.33), (0.57, 0.33), color=RED, dashed=True)
    arrow((0.87, 0.58), (0.80, 0.45), color=GREY)
    # Keep the interpretation rule inside a dedicated footer container rather
    # than as free-floating text below the diagram.  This prevents apparent
    # overflow when the figure is scaled by Word or a submission portal.
    box(0.03, 0.015, 0.94, 0.18, "Claim ceiling",
        ["Wet-lab handoff: identity -> preparation -> exposure -> intervention",
         "A high rank prioritizes a test; it does not establish a biological mechanism."],
        RED)
    ax.text(0.005, 0.98, "A", fontsize=12, fontweight="bold", va="top")
    _save(fig, "Figure1_framework_and_temporal_design")


def figure2_coverage():
    protocol = json.loads((BENCH / "protocol.json").read_text(encoding="utf-8"))
    manifest = json.loads((BENCH / "candidate_corpus.json").read_text(encoding="utf-8"))
    chembl = json.loads((BENCH / "cutoff_chembl_2022.json").read_text(encoding="utf-8"))
    formulas = protocol["formulas"]
    train = [len(manifest["formulas"][f]["train"]["pmids"]) for f in formulas]
    test = [len(manifest["formulas"][f]["test"]["pmids"]) for f in formulas]
    coverage = []
    labels = []
    for f in formulas:
        compounds = {c for cs in chembl["formula_compounds"][f].values() for c in cs}
        available = sum(chembl["compounds"][c]["status"] == "available"
                        for c in compounds)
        coverage.append(100 * available / len(compounds))
        labels.append(f"{available}/{len(compounds)}")

    order = np.argsort(np.array(test) + np.array(train))
    y = np.arange(len(formulas))
    names = [manifest["formulas"][formulas[i]]["aliases"][0] for i in order]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.1, 4.85),
                                   gridspec_kw={"width_ratios": [1.25, 1]})
    train_o = np.array(train)[order]
    test_o = np.array(test)[order]
    ax1.barh(y, train_o, color=BLUE, label="Pre-cutoff (<=2022)")
    ax1.barh(y, test_o, left=train_o, color=ORANGE,
             label="Held-out screening ledger (2023-2025)")
    ax1.set_yticks(y, names)
    ax1.set_xlabel("Retrieved PubMed articles, n")
    ax1.set_title("A  Temporal corpus (n = 932 records)", pad=8)
    ax1.legend(frameon=False, loc="lower right", fontsize=7)
    ax1.grid(axis="x", color=LIGHT, linewidth=0.7)
    ax1.set_axisbelow(True)
    totals = train_o + test_o
    for total, yy in zip(totals, y):
        ax1.text(total + 1.2, yy, f"{total}", va="center", fontsize=6.7,
                 color=DARK, ha="left")
    for s in ("top", "right"):
        ax1.spines[s].set_visible(False)

    cov_o = np.array(coverage)[order]
    lab_o = np.array(labels)[order]
    ax2.hlines(y, 0, cov_o, color=LIGHT, linewidth=4)
    ax2.scatter(cov_o, y, color=TEAL, s=28, zorder=3)
    for value, yy, label in zip(cov_o, y, lab_o):
        ax2.text(min(value + 2, 94), yy, label, va="center", fontsize=6.7,
                 color=DARK, ha="left" if value < 92 else "right")
    ax2.set_yticks(y, [])
    ax2.set_xlim(0, 100)
    ax2.set_xlabel("Compounds with eligible ChEMBL activity, %")
    ax2.set_title("B  Measured coverage (37/67 eligible)", pad=8)
    ax2.grid(axis="x", color=LIGHT, linewidth=0.7)
    ax2.set_axisbelow(True)
    for s in ("top", "right", "left"):
        ax2.spines[s].set_visible(False)
    fig.text(0.01, 0.01,
             "ChEMBL: current database, human targets, document year <=2022; no historical snapshot claim.",
             fontsize=6.7, color=GREY)
    fig.tight_layout(rect=(0, 0.035, 1, 1))
    _save(fig, "Figure2_corpus_and_measured_coverage")


def figure_s1_screening_ledger():
    """Show the auditable screening workflow without implying performance."""
    corpus = json.loads((BENCH / "candidate_corpus.json").read_text(encoding="utf-8"))
    formulas = json.loads((BENCH / "protocol.json").read_text(encoding="utf-8"))["formulas"]
    categories = ["proxy eligible", "exact/no assay signal", "manual identity review", "machine exclusion"]
    colors = [TEAL, BLUE, ORANGE, RED]
    values = {key: [] for key in categories}
    labels = []
    for formula in formulas:
        test = corpus["formulas"][formula]["test"]
        mentions = test["mentions"]
        eligible = set(test["proxy_eligible_pmids"])
        values["proxy eligible"].append(len(eligible))
        values["exact/no assay signal"].append(sum(
            item["eligibility"]["status"] == "eligible"
            and item["triage"]["status"] != "candidate_experimental"
            for item in mentions.values()))
        values["manual identity review"].append(sum(
            item["eligibility"]["status"] == "manual_review"
            for item in mentions.values()))
        values["machine exclusion"].append(sum(
            item["eligibility"]["status"] == "exclude"
            for item in mentions.values()))
        labels.append(corpus["formulas"][formula]["aliases"][0])
    order = np.argsort(np.array([sum(values[key][i] for key in categories)
                                 for i in range(len(formulas))]))
    y = np.arange(len(formulas))
    fig, ax = plt.subplots(figsize=(7.1, 4.7))
    left = np.zeros(len(formulas))
    for category, color in zip(categories, colors):
        series = np.array(values[category])[order]
        ax.barh(y, series, left=left, label=category, color=color)
        left += series
    ax.set_yticks(y, np.array(labels)[order])
    ax.set_xlabel("Machine-retrieved 2023-2025 records, n")
    ax.set_title("Held-out screening ledger (n = 339)", pad=8)
    ax.legend(frameon=False, fontsize=7, ncol=2, loc="lower right")
    ax.grid(axis="x", color=LIGHT, linewidth=0.7)
    ax.set_axisbelow(True)
    totals = left.copy()
    for total, yy in zip(totals, y):
        ax.text(total + 0.7, yy, f"{int(total)}", va="center", fontsize=6.6,
                color=DARK, ha="left")
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    fig.text(0.01, 0.01,
             "Eligibility is a screening rule, not a reference standard; all final inclusion requires blinded full-text review.",
             fontsize=6.8, color=RED, fontweight="bold")
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    _save(fig, "FigureS1_screening_ledger")


def _evidence_composition(feature: dict) -> str:
    """Classify the evidence sources visible for one ranked target."""
    measured = bool(feature["strong_compounds"] or feature["weak_compounds"])
    prediction = bool(feature["stp_compounds"])
    literature = bool(feature["literature_pmids"])
    n = int(measured) + int(prediction) + int(literature)
    if n >= 2:
        return "mixed"
    if measured:
        return "measured-only"
    if prediction:
        return "prediction-only"
    if literature:
        return "literature-only"
    return "no recorded evidence"


def figure3_counterfactual_audit():
    """Audit how evidence rules change decisions without claiming accuracy."""
    protocol = json.loads((BENCH / "protocol.json").read_text(encoding="utf-8"))
    corpus = json.loads((BENCH / "candidate_corpus.json").read_text(encoding="utf-8"))
    rankings = json.loads((BENCH / "locked_rankings_and_proxy_audit.json").read_text(encoding="utf-8"))
    chembl = json.loads((BENCH / "cutoff_chembl_2022.json").read_text(encoding="utf-8"))

    rows = []
    composition = {"ungated_union": Counter(), "evidence_gated": Counter()}
    for formula in protocol["formulas"]:
        entry = rankings["formulas"][formula]
        gated = entry["rankings"]["evidence_gated"]
        for baseline, label in (("ungated_union", "ungated_union"),
                                ("cutoff_chembl_measured_only", "measured_only")):
            base = entry["rankings"][baseline]
            k10 = min(10, len(gated), len(base))
            gated10, base10 = set(gated[:k10]), set(base[:k10])
            k20 = min(20, len(gated), len(base))
            gated20, base20 = set(gated[:k20]), set(base[:k20])
            common = gated20 & base20
            displacement = [abs((gated.index(target) + 1) -
                                (base.index(target) + 1)) for target in common]
            rows.append({
                "formula": formula,
                "label": corpus["formulas"][formula]["aliases"][0],
                "baseline": label,
                "top_k": k10,
                "top10_jaccard": (len(gated10 & base10) / len(gated10 | base10)
                                   if gated10 | base10 else 1.0),
                "top20_candidate_loss": len(base20 - gated20),
                "top20_common": len(common),
                "mean_abs_rank_displacement": (sum(displacement) / len(displacement)
                                                if displacement else 0.0),
            })
        for method in ("ungated_union", "evidence_gated"):
            for target in entry["rankings"][method][:20]:
                composition[method][_evidence_composition(entry["features"][target])] += 1

    def summarize(label):
        values = [r for r in rows if r["baseline"] == label]
        j = [r["top10_jaccard"] for r in values]
        loss = [r["top20_candidate_loss"] for r in values]
        disp = [r["mean_abs_rank_displacement"] for r in values]
        return {
            "median_top10_jaccard": float(np.median(j)),
            "mean_top10_jaccard": float(np.mean(j)),
            "minimum_top10_jaccard": float(np.min(j)),
            "mean_top20_candidate_loss": float(np.mean(loss)),
            "median_top20_candidate_loss": float(np.median(loss)),
            "mean_abs_rank_displacement": float(np.mean(disp)),
        }

    summary = {
        "ungated_union_vs_evidence_gated": summarize("ungated_union"),
        "measured_only_vs_evidence_gated": summarize("measured_only"),
        "top20_slots_per_method": len(protocol["formulas"]) * 20,
        "composition_counts": {k: dict(v) for k, v in composition.items()},
    }

    # Panel A: paired overlap across the 14 formulas.
    u_rows = [r for r in rows if r["baseline"] == "ungated_union"]
    m_rows = [r for r in rows if r["baseline"] == "measured_only"]
    order = np.argsort([r["top10_jaccard"] for r in u_rows])
    names = [u_rows[i]["label"] for i in order]
    y = np.arange(len(names))
    fig, (ax1, ax2, ax3) = plt.subplots(
        1, 3, figsize=(7.1, 4.35),
        gridspec_kw={"width_ratios": [1.35, 1.0, 1.15]},
    )
    u_j = np.array([u_rows[i]["top10_jaccard"] for i in order])
    m_by_formula = {r["formula"]: r for r in m_rows}
    m_j = np.array([m_by_formula[u_rows[i]["formula"]]["top10_jaccard"] for i in order])
    ax1.plot(u_j, y, "o", color=BLUE, markersize=4, label="Ungated union")
    ax1.plot(m_j, y, "D", color=ORANGE, markersize=3.5, label="Measured-only")
    ax1.set_yticks(y, names)
    ax1.tick_params(axis="y", labelsize=5.9)
    ax1.set_xlim(0, 1.05)
    ax1.set_xlabel("Top-10 Jaccard with gated rank")
    ax1.set_title("A  Ranking overlap", pad=8)
    ax1.grid(axis="x", color=LIGHT, linewidth=0.7)
    ax1.set_axisbelow(True)
    ax1.legend(frameon=False, fontsize=5.9, loc="upper center",
              bbox_to_anchor=(0.5, -0.16), ncol=2)
    for s in ("top", "right"):
        ax1.spines[s].set_visible(False)

    categories = ["measured-only", "mixed", "literature-only"]
    colors = [BLUE, TEAL, ORANGE]
    methods = [("ungated_union", "Ungated union"), ("evidence_gated", "Evidence-gated")]
    left = np.zeros(2)
    for category, color in zip(categories, colors):
        vals = np.array([composition[m].get(category, 0) for m, _ in methods])
        ax2.barh(np.arange(2), vals, left=left, color=color, label=category)
        for i, (start, value) in enumerate(zip(left, vals)):
            if value >= 12:
                ax2.text(start + value / 2, i, f"{int(value)}", ha="center",
                         va="center", fontsize=6.3, color="white", fontweight="bold")
        left += vals
    ax2.set_yticks(np.arange(2), [label for _, label in methods])
    ax2.tick_params(axis="y", labelsize=6.1)
    ax2.set_xlim(0, 280)
    ax2.set_xlabel("Ranked target slots, n")
    ax2.set_title("B  Top-20 composition (n = 280)", pad=8)
    ax2.legend(frameon=False, fontsize=5.6, loc="upper center",
               bbox_to_anchor=(0.5, -0.16), ncol=3)
    ax2.grid(axis="x", color=LIGHT, linewidth=0.7)
    ax2.set_axisbelow(True)
    for s in ("top", "right"):
        ax2.spines[s].set_visible(False)
    ax2.text(0, -0.34, "Prediction-only = 0 in the fixed primary universe",
             transform=ax2.transAxes, fontsize=5.8, color=GREY)

    # Panel C: observed compound-level attrition, followed by unmeasured claim gates.
    statuses = chembl["compounds"]
    total = len(statuses)
    resolved = sum(bool(v.get("chembl_id")) for v in statuses.values())
    eligible = sum(v.get("status") == "available" for v in statuses.values())
    vals = [total, resolved, eligible]
    y3 = np.arange(3)
    ax3.barh(y3, vals, color=[ORANGE, BLUE, TEAL], height=0.55)
    ax3.set_yticks(y3, ["Marker compounds", "Resolved IDs", "Temporal measured activity"])
    ax3.tick_params(axis="y", labelsize=6.0)
    ax3.invert_yaxis()
    ax3.set_xlim(0, 75)
    ax3.set_xlabel("Compounds, n")
    ax3.set_title("C  Attrition and claim gates", pad=8)
    for yy, value, pct in zip(y3, vals, [100.0, 100 * resolved / total, 100 * eligible / total]):
        ax3.text(value + 1.2, yy, f"{value} ({pct:.1f}%)", va="center", fontsize=6.5,
                 color=DARK)
    ax3.text(0, -0.32,
             "Preparation detection → exposure → selective intervention\n"
             "not assessed in this study (not zero)",
             transform=ax3.transAxes, fontsize=5.8, color=RED)
    ax3.grid(axis="x", color=LIGHT, linewidth=0.7)
    ax3.set_axisbelow(True)
    for s in ("top", "right"):
        ax3.spines[s].set_visible(False)

    fig.text(0.01, 0.015,
             "Counterfactual rankings use the same fixed candidate universe; overlap and displacement describe decision sensitivity, not biological accuracy.",
             fontsize=6.2, color=GREY)
    fig.tight_layout(rect=(0, 0.16, 1, 1))
    _save(fig, "Figure3_counterfactual_evidence_audit")

    source = {
        "protocol_id": protocol["protocol_id"],
        "definitions": {
            "ungated_union": "existing weighted union rank without the lexical evidence-tier gate",
            "evidence_gated": "existing deterministic rank with the evidence-tier gate",
            "measured_only": "pre-cutoff ChEMBL measured-only rank; no literature or STP tie-break",
            "top10_jaccard": "set intersection divided by set union for the available top-10 lists",
            "top20_candidate_loss": "baseline top-20 targets absent from evidence-gated top-20",
            "composition": "presence/absence of measured, STP, and pre-cutoff literature fields",
        },
        "per_formula": rows,
        "summary": summary,
        "compound_attrition": {
            "marker_compounds": total,
            "resolved_chembl_ids": resolved,
            "temporal_measured_eligible": eligible,
            "preparation_detection": "not_assessed",
            "exposure": "not_assessed",
            "selective_intervention": "not_assessed",
        },
    }
    (OUT / "figure3_counterfactual_audit.json").write_text(
        json.dumps(source, ensure_ascii=False, indent=2), encoding="utf-8")


def figure4_sensitivity():
    """Summarise threshold and source-ablation sensitivity without accuracy claims."""
    def read_tsv(name):
        with (AUDITS / name).open(encoding="utf-8") as f:
            return list(csv.DictReader(f, delimiter="\t"))
    threshold = read_tsv("threshold_sensitivity.tsv")
    turnover = read_tsv("top10_turnover.tsv")
    ablation = read_tsv("source_ablation.tsv")
    levels = sorted({float(r["pchembl_threshold"]) for r in threshold})
    med_j = [float(np.median([float(r["top10_jaccard_vs_primary"]) for r in threshold if float(r["pchembl_threshold"]) == x])) for x in levels]
    methods = ["evidence_gated_vs_ungated_union", "evidence_gated_vs_cutoff_chembl_measured_only"]
    labels = ["Ungated union", "Measured-only"]
    med_repl = [float(np.median([int(r["top10_replacements"]) for r in turnover if r["comparison"] == x])) for x in methods]
    ab_names = ["primary", "without_literature", "without_measured", "without_measured_and_literature", "literature_only", "measured_only"]
    ab_labels = ["Primary", "− literature", "− measured", "− both", "Literature only", "Measured only"]
    ab_j = [float(np.median([float(r["top10_jaccard_vs_primary"]) for r in ablation if r["ablation"] == x])) for x in ab_names]
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(7.1, 3.0), gridspec_kw={"width_ratios": [1.15, 0.95, 1.25]})
    ax1.plot(levels, med_j, marker="o", color=BLUE, linewidth=1.7, markersize=3.5)
    ax1.axvline(5, color=RED, linestyle="--", linewidth=0.9)
    ax1.set_xlabel("pChEMBL threshold")
    ax1.set_ylabel("Median Top-10 Jaccard vs primary")
    ax1.set_title("A  Threshold sensitivity", pad=8)
    ax1.set_ylim(0, 1.05); ax1.grid(axis="y", color=LIGHT, linewidth=0.7); ax1.set_axisbelow(True)
    for s in ("top", "right"): ax1.spines[s].set_visible(False)
    ax2.bar(np.arange(2), med_repl, color=[TEAL, ORANGE], width=0.55)
    ax2.set_xticks(np.arange(2), labels, rotation=25, ha="right", fontsize=6.2)
    ax2.set_ylabel("Median Top-10 replacements")
    ax2.set_title("B  Decision turnover", pad=8)
    ax2.set_ylim(0, max(med_repl + [1]) + 1); ax2.grid(axis="y", color=LIGHT, linewidth=0.7); ax2.set_axisbelow(True)
    for x, v in enumerate(med_repl): ax2.text(x, v + 0.12, f"{v:.1f}", ha="center", fontsize=7)
    for s in ("top", "right"): ax2.spines[s].set_visible(False)
    ypos = np.arange(len(ab_labels))
    ax3.barh(ypos, ab_j, color=[BLUE, ORANGE, TEAL, RED, ORANGE, BLUE], alpha=0.9)
    ax3.set_yticks(ypos, ab_labels, fontsize=6.1); ax3.invert_yaxis(); ax3.set_xlim(0, 1.05)
    ax3.set_xlabel("Median Top-10 Jaccard vs primary")
    ax3.set_title("C  Source ablation", pad=8)
    ax3.grid(axis="x", color=LIGHT, linewidth=0.7); ax3.set_axisbelow(True)
    for s in ("top", "right"): ax3.spines[s].set_visible(False)
    fig.text(0.01, 0.01, "Sensitivity and ablation quantify ranking dependence; they do not estimate predictive performance.", fontsize=6.2, color=GREY)
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    _save(fig, "Figure4_sensitivity_and_source_ablation")


def main():
    figure1_framework()
    figure2_coverage()
    figure_s1_screening_ledger()
    figure3_counterfactual_audit()
    figure4_sensitivity()
    hashes = {}
    for path in sorted(OUT.iterdir()):
        if path.is_file() and path.name != "figure_manifest.json":
            hashes[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
    manifest = {
        "generated_from": [
            "benchmarks/papera_temporal_v2/protocol.json",
            "benchmarks/papera_temporal_v2/candidate_corpus.json",
            "benchmarks/papera_temporal_v2/cutoff_chembl_2022.json",
            "benchmarks/papera_temporal_v2/locked_rankings_and_proxy_audit.json",
        ],
        "files_sha256": hashes,
    }
    (OUT / "figure_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"[OK] figures -> {OUT}")


if __name__ == "__main__":
    main()
