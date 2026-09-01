"""ChEMBL REST API client (read-only, no key required)."""
from __future__ import annotations

import json
import time
from typing import List, Optional

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

DEFAULT_BASE = "https://www.ebi.ac.uk/chembl/api/data"

# 批量查询时 ChEMBL 会限流（429），用指数退避自动重试，
# 避免一次复现任务因瞬时限流产生大量 api_failure。
_RETRY = Retry(
    total=4,
    backoff_factor=1.0,
    status_forcelist=[429, 500, 502, 503, 504],
    allowed_methods=["GET"],
    respect_retry_after_header=True,
)


def _as_float(value) -> Optional[float]:
    """Coerce ChEMBL numeric strings to float; None if not numeric."""
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


class ChEMBLClient:
    def __init__(self, base_url: str = DEFAULT_BASE, timeout_s: int = 20,
                 offline: bool = False):
        self.base = base_url.rstrip("/")
        self.timeout = timeout_s
        self.offline = offline
        self._session = requests.Session()
        adapter = HTTPAdapter(max_retries=_RETRY)
        self._session.mount("https://", adapter)
        self._session.mount("http://", adapter)

    def search_molecules(self, query: str, limit: int = 5) -> List[dict]:
        """Search molecules by name/synonym. Returns list of molecule dicts."""
        if self.offline:
            return []
        url = f"{self.base}/molecule/search"
        params = {"q": query, "format": "json", "limit": limit}
        resp = self._session.get(url, params=params, timeout=self.timeout)
        resp.raise_for_status()
        return resp.json().get("molecules", [])

    def target_activities(self, chembl_id: str, limit: int = 50,
                          organism: Optional[str] = None,
                          order_by: Optional[str] = None,
                          document_year_lte: Optional[int] = None,
                          document_year_gte: Optional[int] = None) -> List[dict]:
        """Return bioactivities for a molecule against targets.

        Optional filters (behavior unchanged when omitted):
        - organism: restrict to one target organism, e.g. "Homo sapiens"
        - order_by: ChEMBL field, e.g. "-pchembl_value" for potency-descending
        - document_year_lte/gte: restrict records by the associated document
          year. These server-side filters are used by temporal benchmarks to
          prevent post-cutoff activity records from leaking into training.
        """
        if self.offline:
            return []
        url = f"{self.base}/activity"
        params = {"molecule_chembl_id": chembl_id, "format": "json",
                  "limit": limit, "pchembl_value__isnull": "false"}
        if organism:
            params["target_organism"] = organism
        if order_by:
            params["order_by"] = order_by
        if document_year_lte is not None:
            params["document_year__lte"] = int(document_year_lte)
        if document_year_gte is not None:
            params["document_year__gte"] = int(document_year_gte)
        resp = self._session.get(url, params=params, timeout=self.timeout)
        resp.raise_for_status()
        activities = resp.json().get("activities", [])
        out = []
        for act in activities:
            tgt = act.get("target_chembl_id")
            if not tgt:
                continue
            out.append({
                "target_chembl_id": tgt,
                "target_pref_name": act.get("target_pref_name"),
                "standard_type": act.get("standard_type"),
                "standard_value": act.get("standard_value"),
                "standard_units": act.get("standard_units"),
                "pchembl_value": _as_float(act.get("pchembl_value")),
                "document_chembl_id": act.get("document_chembl_id"),
                "document_year": act.get("document_year"),
                "assay_chembl_id": act.get("assay_chembl_id"),
                "assay_type": act.get("assay_type"),
            })
        return out

    def resolve_compound_targets(self, compound_names: List[str],
                                 sleep_s: float = 0.35,
                                 organism: Optional[str] = None,
                                 min_pchembl_value: Optional[float] = None,
                                 order_by: Optional[str] = None,
                                 document_year_lte: Optional[int] = None,
                                 document_year_gte: Optional[int] = None,
                                 top_n: int = 30,
                                 fetch_limit: int = 100) -> dict:
        """Map compound names -> ChEMBL IDs -> target preferred names.

        Activities are potency-ranked (pchembl_value desc, or order_by if
        given) and truncated to top_n; min_pchembl_value keeps only records
        with pchembl_value (log potency) >= the threshold, e.g. 5.0 (~10 uM).
        Deterministic ordering; failures are recorded, never fabricated.
        """
        result: dict = {}
        for name in compound_names:
            try:
                mols = self.search_molecules(name, limit=1)
                if not mols:
                    result[name] = {"chembl_id": None, "targets": []}
                    continue
                chembl_id = mols[0]["molecule_chembl_id"]
                acts = self.target_activities(chembl_id, limit=fetch_limit,
                                              organism=organism, order_by=order_by,
                                              document_year_lte=document_year_lte,
                                              document_year_gte=document_year_gte)
                if min_pchembl_value is not None:
                    acts = [a for a in acts
                            if (a.get("pchembl_value") or 0.0) >= min_pchembl_value]
                acts = sorted(acts, key=lambda a: a.get("pchembl_value") or 0.0,
                              reverse=True)[:top_n]
                targets = sorted({a["target_pref_name"] for a in acts
                                  if a.get("target_pref_name")})
                result[name] = {"chembl_id": chembl_id, "targets": targets}
                time.sleep(sleep_s)
            except (requests.RequestException, ValueError):
                result[name] = {"chembl_id": None, "targets": [], "error": "api_failure"}
        return result
