import importlib.util
import json
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "make_papera_figures.py"
SPEC = importlib.util.spec_from_file_location("make_papera_figures", SCRIPT)
FIGURES = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(FIGURES)


def test_top10_replacement_summary_uses_only_complete_top10_lists():
    rankings = json.loads(
        (FIGURES.BENCH / "locked_rankings_and_proxy_audit.json").read_text(
            encoding="utf-8"
        )
    )

    result = FIGURES.top10_replacement_summary(
        rankings["formulas"], list(rankings["formulas"])
    )

    assert result["evidence_gated_vs_ungated_union"] == {
        "n": 14,
        "median_replacements": 1.0,
    }
    assert result["evidence_gated_vs_cutoff_chembl_measured_only"] == {
        "n": 13,
        "median_replacements": 4.0,
    }


def test_short_measured_only_list_is_reported_at_available_k_not_as_top10():
    rankings = json.loads(
        (FIGURES.BENCH / "locked_rankings_and_proxy_audit.json").read_text(
            encoding="utf-8"
        )
    )

    rows = FIGURES.available_k_turnover_rows(
        rankings["formulas"], list(rankings["formulas"])
    )
    short = next(
        row
        for row in rows
        if row["formula"] == "六味地黄丸"
        and row["comparison"] == "evidence_gated_vs_cutoff_chembl_measured_only"
    )

    assert short["available_k"] == 3
    assert short["replacements_at_k"] == 0
    assert "top10_replacements" not in short
    assert len(rows) == 28


def test_figure_generator_labels_no_human_adjudication_as_completed():
    source = SCRIPT.read_text(encoding="utf-8")
    assert "No human-adjudicated labels" in source
    assert "Blinded review + adjudication" not in source
    assert "human full-text adjudication was not performed" in source


def test_generated_svg_has_no_trailing_whitespace():
    FIGURES.figure3_counterfactual_audit()

    svg = (FIGURES.OUT / "Figure3_counterfactual_evidence_audit.svg").read_text(
        encoding="utf-8"
    )
    assert all(line == line.rstrip() for line in svg.splitlines())
