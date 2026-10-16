"""Render the exploratory structure-matched transfer-audit summary.

This script uses only the derived, public summaries in ``transfer_audit``.
It intentionally does not retrieve external records or reconstruct a biological
reference standard.  The output is descriptive workflow-transfer evidence.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
SUMMARY = ROOT / "transfer_audit" / "results_summary.json"
OUT = ROOT / "figures" / "Figure5_transfer_audit.png"


def build_figure(data: dict):
    fig, (ax_a, ax_b, ax_c) = plt.subplots(
        1, 3, figsize=(7.5, 3.8), dpi=400, facecolor="white",
        gridspec_kw={"width_ratios": [0.9, 1.0, 2.1], "wspace": 0.9},
    )
    blue, orange, green, red = "#2f6690", "#d97706", "#3f8f6b", "#b84a62"

    # Panel A: a compact vertical flow keeps each denominator legible.
    ax_a.set_axis_off()
    ax_a.text(0.0, 1.03, "A  Task status", transform=ax_a.transAxes,
              fontsize=8.5, fontweight="bold", color="#243447")
    n = data["n_tasks_selected"]
    r = data["n_tasks_rankable"]
    nr = data["n_tasks_not_rankable"]
    ax_a.text(0.02, 0.87, "All selected tasks completed retrieval", transform=ax_a.transAxes,
              fontsize=7.1, color="#52606d")
    flow = [(0.68, n, "selected", blue), (0.42, r, "rankable", green),
            (0.16, nr, "without eligible candidate", red)]
    for y, value, label, color in flow:
        ax_a.text(0.13, y, f"{value}", transform=ax_a.transAxes, ha="center", va="center",
                  fontsize=13, fontweight="bold", color=color)
        ax_a.text(0.32, y, label, transform=ax_a.transAxes, ha="left", va="center",
                  fontsize=7.3, color="#34495e")
    for y0, y1 in [(0.63, 0.49), (0.37, 0.23)]:
        ax_a.annotate("", xy=(0.13, y1), xytext=(0.13, y0), xycoords="axes fraction",
                      arrowprops={"arrowstyle": "-|>", "lw": 1.25, "color": "#829ab1"})

    ax_b.text(-0.12, 1.03, "B  Top-k changes", transform=ax_b.transAxes,
              fontsize=8.5, fontweight="bold", color="#243447")
    comps = [data["top10_change_vs_ungated"], data["top10_change_vs_measured_only"]]
    labels = ["ungated", "measured\npriority"]
    vals = [c["tasks_with_any_top10_replacement"] for c in comps]
    pcts = [c["percent_tasks_with_replacement"] for c in comps]
    ax_b.bar(labels, vals, color=[blue, orange], width=0.58)
    ax_b.set_ylim(0, max(vals) + 3)
    ax_b.set_ylabel("Tasks with ≥1 replacement", fontsize=7.0, labelpad=4)
    ax_b.tick_params(axis="x", labelsize=6.8, pad=5)
    ax_b.tick_params(axis="y", labelsize=6.8)
    for spine in ax_b.spines.values():
        spine.set_color("#aab7c4")
        spine.set_linewidth(0.8)
    ax_b.grid(axis="y", color="#e5e7eb", linewidth=0.7)
    ax_b.set_axisbelow(True)
    for i, (v, p) in enumerate(zip(vals, pcts)):
        ax_b.text(i, v + 0.35, f"{v}/36 ({p:.1f}%)", ha="center", fontsize=7.1,
                  color="#243447")
    ax_b.text(0.5, -0.29, "k = min(10, N); N = task candidate-universe size",
              transform=ax_b.transAxes, ha="center", fontsize=6.4, color="#52606d",
              wrap=True)
    ax_b.text(0.5, -0.39, "Measured-priority uses the shared candidate universe",
              transform=ax_b.transAxes, ha="center", fontsize=6.4, color="#52606d",
              wrap=True)

    ax_c.text(0.0, 1.03, "C  Evidence-object counts", transform=ax_c.transAxes,
              fontsize=8.5, fontweight="bold", color="#243447")
    measured = data["evidence_object_counts_by_layer_and_status"]
    layers = [
        ("Activity records", data["activity_observation_rows"], blue),
        ("Eligible measured", measured["measured|eligible"], green),
        ("PubMed records", data["pubmed_records_across_tasks_not_independent"], orange),
        ("Screening only", measured["literature|screening_only"], blue),
        ("Manual review", measured["literature|manual_review"], orange),
        ("Excluded", measured["literature|excluded"], red),
    ]
    ax_c.set_xlim(0, 1.38)
    ax_c.set_ylim(-0.65, len(layers) - 0.45)
    max_value = max(value for _, value, _ in layers)
    y_positions = list(range(len(layers)))
    for y, (label, value, color) in zip(y_positions, layers):
        width = value / max_value
        ax_c.barh(y, width, height=0.58, color=color, alpha=0.9)
        ax_c.text(width + 0.025, y, f"{value:,}", va="center", fontsize=7.1,
                  color="#243447")
    ax_c.set_yticks(y_positions, [label for label, _, _ in layers], fontsize=6.5)
    ax_c.invert_yaxis()
    ax_c.set_xticks([])
    ax_c.tick_params(axis="y", length=0, pad=3)
    for spine in ax_c.spines.values():
        spine.set_visible(False)
    ax_c.text(0.0, -0.17,
              "Counts are nested evidence objects; PubMed records can recur across tasks and are not independent.",
              transform=ax_c.transAxes, ha="left", fontsize=6.4, color="#52606d")

    fig.suptitle("Exploratory structure-matched transfer audit", fontsize=10,
                 fontweight="bold", color="#243447", y=0.995)
    fig.subplots_adjust(top=0.84, bottom=0.23, left=0.04, right=0.98)

    return fig


def main() -> None:
    data = json.loads(SUMMARY.read_text(encoding="utf-8"))
    fig = build_figure(data)
    fig.savefig(OUT, bbox_inches="tight", facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
