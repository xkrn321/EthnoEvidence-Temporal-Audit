"""Regression tests for transfer-audit figure semantics."""
import json
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "make_transfer_audit_figure.py"
SPEC = spec_from_file_location("transfer_audit_figure", SCRIPT)
MODULE = module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(MODULE)


def test_transfer_figure_labels_rankable_task_topk_and_shared_measured_ordering():
    summary = json.loads(MODULE.SUMMARY.read_text(encoding="utf-8"))
    fig = MODULE.build_figure(summary)
    try:
        panel_b = fig.axes[1]
        visible_text = [text.get_text() for text in panel_b.texts]
        x_labels = [tick.get_text() for tick in panel_b.get_xticklabels()]
        bar_labels = [text.get_text() for text in panel_b.texts if "/36" in text.get_text()]

        assert "B  Top-k changes" in visible_text
        assert "k = min(10, N); N = task candidate-universe size" in visible_text
        assert x_labels == ["ungated", "measured\npriority"]
        assert bar_labels == ["6/36 (16.7%)", "9/36 (25.0%)"]
    finally:
        plt.close(fig)
