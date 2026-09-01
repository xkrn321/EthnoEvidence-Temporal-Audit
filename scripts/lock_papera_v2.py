"""Write a reproducibility checksum manifest for the v2 release inputs.

The manifest is deliberately limited to the protocol, raw PubMed records,
screening ledger, pharmacology tier, and locked ranking output that support the
manuscript.  It is not a DOI and does not turn mutable upstream services into a
historical snapshot.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BENCH = ROOT / "benchmarks" / "papera_temporal_v2"
OUTPUT = BENCH / "LOCK.sha256.json"
CORE = (
    "protocol.json",
    "PROTOCOL.md",
    "RANKING_SPEC.md",
    "DEVIATIONS.md",
    "candidate_corpus.json",
    "manual_screening.tsv",
    "cutoff_chembl_2022.json",
    "locked_rankings_and_proxy_audit.json",
    "retrieval_manifest.json",
)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    paths = [BENCH / name for name in CORE]
    paths.extend(sorted((BENCH / "raw_pubmed").glob("*_all.json")))
    missing = [path for path in paths if not path.is_file()]
    if missing:
        raise SystemExit("Missing release inputs: " + ", ".join(map(str, missing)))

    release = {
        "release_id": "papera_temporal_v2",
        "checksum_algorithm": "sha256",
        "files": [
            {"path": str(path.relative_to(BENCH)), "sha256": digest(path)}
            for path in paths
        ],
    }
    release["file_count"] = len(release["files"])
    OUTPUT.write_text(json.dumps(release, ensure_ascii=False, indent=2) + "\n",
                      encoding="utf-8")
    print(OUTPUT)


if __name__ == "__main__":
    main()
