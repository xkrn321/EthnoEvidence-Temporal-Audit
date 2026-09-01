"""Regression tests for Paper A temporal-ranking controls."""
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "build_temporal_rankings.py"
SPEC = spec_from_file_location("papera_rankings", SCRIPT)
MODULE = module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(MODULE)


def _feature(strong=(), weak=(), stp=(), lit=(), pchembl=0.0, probability=0.0):
    return {
        "strong_compounds": set(strong),
        "weak_compounds": set(weak),
        "stp_compounds": set(stp),
        "literature_pmids": set(lit),
        "best_pchembl": pchembl,
        "best_stp_probability": probability,
    }


def test_measured_only_cannot_use_stp_or_literature_tiebreaks():
    features = {
        "ZGENE": _feature(weak=["a"], pchembl=4.5, stp=["a"], lit=["1", "2"]),
        "AGENE": _feature(weak=["b"], pchembl=4.5),
    }
    assert MODULE._measured_only(features) == ["AGENE", "ZGENE"]


def test_primary_rank_ignores_current_stp_signal():
    features = {
        "GENEA": _feature(lit=["1"], stp=["a"], probability=0.99),
        "GENEB": _feature(weak=["b"], pchembl=4.2),
    }
    assert MODULE._sorted(features, include_stp=False)[0] == "GENEB"
