import importlib.util
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "verify_papera_v2_release.py"
SPEC = importlib.util.spec_from_file_location("verify_papera_v2_release", SCRIPT)
VERIFIER = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(VERIFIER)


def test_verifier_describes_future_records_as_screening_not_holdout(capsys):
    assert VERIFIER.main() == 0
    output = capsys.readouterr().out

    assert "post-cutoff screening (2023-2025)" in output
    assert "held-out (2023-2025)" not in output
