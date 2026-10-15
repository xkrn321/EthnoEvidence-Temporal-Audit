"""Render the exploratory structure-matched transfer-audit summary.

This script uses only the derived, public summaries in ``transfer_audit``.
It intentionally does not retrieve external records or reconstruct a biological
reference standard.  The output is descriptive workflow-transfer evidence.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch


ROOT = Path(__file__).resolve().parents[1]
SUMMARY = ROOT / "transfer_audit" / "results_summary.json"
OUT = ROOT / "figures" / "Figure5_transfer_audit.png"


def box(ax, x, y, w, h, label, value, color):
    patch = FancyBboxPatch(
        (x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.02",
        linewidth=1.3, edgecolor=color, facecolor="#f7fafc",
    )
    ax.add_patch(patch)
    ax.text(x + w / 2, y + h * 0.62, value, ha="center", va="center",
            fontsize=18, fontweight="bold", color="#243447")
    ax.text(x + w / 2, y + h * 0.28, label, ha="center", va="center",
            fontsize=9.5, color="#52606d", wrap=True)


def arrow(ax, x0, x1, y):
    ax.annotate("", xy=(x1, y), xytext=(x0, y),
                arrowprops={"arrowstyle": "-|>", "lw": 1.5, "color": "#829ab1"})


def main() -> None:
    data = json.loads(SUMMARY.read_text(encoding="utf-8"))
    fig = plt.figure(figsize=(11, 5.75), dpi=200, facecolor="white")
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.05], hspace=0.48,
                          wspace=0.26)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, :])
    for ax in (ax_a, ax_b, ax_c):
        ax.set_axis_off()

    blue, orange, green, red = "#2f6690", "#d97706", "#3f8f6b", "#b84a62"
    ax_a.text(0.0, 1.06, "A  Task retrieval and rankability", transform=ax_a.transAxes,
              fontsize=11, fontweight="bold", color="#243447")
    n = data["n_tasks_selected"]
    r = data["n_tasks_rankable"]
    nr = data["n_tasks_not_rankable"]
    box(ax_a, 0.02, 0.32, 0.27, 0.42, "selected tasks", str(n), blue)
    box(ax_a, 0.36, 0.32, 0.27, 0.42, "rankable", f"{r}\n(90.0%)", green)
    box(ax_a, 0.70, 0.32, 0.27, 0.42, "no eligible candidate", str(nr), red)
    arrow(ax_a, 0.29, 0.35, 0.53)
    arrow(ax_a, 0.63, 0.69, 0.53)
    ax_a.text(0.5, 0.08, "All 40 retrievals completed under the frozen eligibility rules.",
              transform=ax_a.transAxes, ha="center", fontsize=8.7, color="#52606d")

    ax_b.text(0.0, 1.06, "B  Top-10 changes among rankable tasks", transform=ax_b.transAxes,
              fontsize=11, fontweight="bold", color="#243447")
    comps = [data["top10_change_vs_ungated"], data["top10_change_vs_measured_only"]]
    labels = ["vs ungated", "vs measured-only"]
    vals = [c["tasks_with_any_top10_replacement"] for c in comps]
    pcts = [c["percent_tasks_with_replacement"] for c in comps]
    ax_b.bar(labels, vals, color=[blue, orange], width=0.56)
    ax_b.set_ylim(0, max(vals) + 4)
    ax_b.set_ylabel("tasks with ≥1 replacement", fontsize=8.5)
    ax_b.tick_params(labelsize=8.5)
    for spine in ax_b.spines.values():
        spine.set_visible(False)
    ax_b.grid(axis="y", color="#e5e7eb", linewidth=0.7)
    ax_b.set_axisbelow(True)
    for i, (v, p) in enumerate(zip(vals, pcts)):
        ax_b.text(i, v + 0.35, f"{v}/36 ({p:.1f}%)", ha="center", fontsize=9,
                  color="#243447")
    ax_b.text(0.5, -0.28, "Median replacements = 0 in both fixed-universe comparisons.",
              transform=ax_b.transAxes, ha="center", fontsize=8.7, color="#52606d")

    ax_c.text(0.0, 1.06, "C  Evidence-object attrition across layers", transform=ax_c.transAxes,
              fontsize=11, fontweight="bold", color="#243447")
    measured = data["evidence_object_counts_by_layer_and_status"]
    layers = [
        ("Activity observations", data["activity_observation_rows"], blue),
        ("Measured eligible", measured["measured|eligible"], green),
        ("PubMed records", data["pubmed_records_across_tasks_not_independent"], orange),
        ("Screening-only", measured["literature|screening_only"], blue),
        ("Manual review", measured["literature|manual_review"], orange),
        ("Excluded", measured["literature|excluded"], red),
    ]
    x_positions = [0.02, 0.18, 0.36, 0.54, 0.70, 0.86]
    widths = [0.13] * 6
    for (label, value, color), x, w in zip(layers, x_positions, widths):
        box(ax_c, x, 0.29, w, 0.48, label, f"{value:,}", color)
    ax_c.text(0.5, 0.04,
              "Counts are nested evidence objects and not independent observations;"
              " licensed raw tables and retrieval caches are not redistributed.",
              transform=ax_c.transAxes, ha="center", fontsize=8.6, color="#52606d")

    fig.suptitle("Exploratory structure-matched transfer audit", fontsize=15,
                 fontweight="bold", color="#243447", y=0.99)
    fig.savefig(OUT, bbox_inches="tight", facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
