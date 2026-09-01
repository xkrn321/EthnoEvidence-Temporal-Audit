"""PubChem PUG REST client (read-only, no key required)."""
from __future__ import annotations

import time
from typing import List, Optional

import requests

DEFAULT_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"


class PubChemClient:
    def __init__(self, base_url: str = DEFAULT_BASE, timeout_s: int = 20,
                 offline: bool = False):
        self.base = base_url.rstrip("/")
        self.timeout = timeout_s
        self.offline = offline

    def name_to_cid(self, name: str) -> Optional[int]:
        """Resolve a compound name to a PubChem CID."""
        if self.offline:
            return None
        url = f"{self.base}/compound/name/{requests.utils.quote(name)}/cids/JSON"
        resp = requests.get(url, timeout=self.timeout)
        if resp.status_code != 200:
            return None
        data = resp.json()
        try:
            return int(data["IdentifierList"]["CID"][0])
        except (KeyError, IndexError, ValueError):
            return None

    def resolve_cids(self, compound_names: List[str]) -> dict:
        """Map compound names -> PubChem CID (deterministic order)."""
        out: dict = {}
        for name in compound_names:
            try:
                out[name] = {"pubchem_cid": self.name_to_cid(name)}
            except requests.RequestException:
                out[name] = {"pubchem_cid": None, "error": "api_failure"}
        return out

    def get_smiles(self, name: str, retries: int = 3) -> Optional[str]:
        """Compound name -> canonical SMILES (via PubChem), retrying on
        server-busy (PUGREST.ServerBusy). Returns None if unresolvable.

        Results are cached on disk (configs/cache/) so repeated analyses
        resolve instantly and never re-hit PubChem.
        """
        if self.offline:
            return None
        from .cache import cache_get, cache_set
        cached = cache_get("smiles", name)
        if cached is not None:
            return cached.get("smiles") if isinstance(cached, dict) else None
        cid = self.name_to_cid(name)
        if not cid:
            cache_set("smiles", name, {"smiles": None})
            return None
        url = (f"{self.base}/compound/cid/{cid}/property/"
               "IsomericSMILES,CanonicalSMILES,ConnectivitySMILES/JSON")
        for attempt in range(retries):
            try:
                r = requests.get(url, timeout=self.timeout)
                if r.status_code != 200:
                    time.sleep(2.0 * (attempt + 1))
                    continue
                props = r.json()["PropertyTable"]["Properties"][0]
                smiles = (props.get("IsomericSMILES")
                          or props.get("CanonicalSMILES")
                          or props.get("ConnectivitySMILES"))
                cache_set("smiles", name, {"smiles": smiles})
                return smiles
            except (requests.RequestException, KeyError, IndexError, ValueError):
                time.sleep(2.0 * (attempt + 1))
        cache_set("smiles", name, {"smiles": None})
        return None

    def similarity_support(self, name: str, threshold: int = 85,
                           max_records: int = 30) -> dict:
        """Structural-support evidence: PubChem 2D Tanimoto similarity search.

        Returns the PubChem CIDs with Tanimoto similarity >= threshold/100.
        Note: PubChem's per-compound gene-target annotation endpoints are not
        exposed via PUG REST (400/Invalid operation), so this modal returns
        structural-support metadata only; it contributes no target claims.
        """
        if self.offline:
            return {"similar_cids": [], "note": "offline"}
        cid = self.name_to_cid(name)
        if not cid:
            return {"similar_cids": [], "note": "no_cid"}
        try:
            r = requests.get(
                f"{self.base}/compound/cid/{cid}/property/IsomericSMILES,CanonicalSMILES,ConnectivitySMILES/JSON",
                timeout=self.timeout)
            if r.status_code != 200:
                return {"similar_cids": [], "note": "no_smiles"}
            props = r.json()["PropertyTable"]["Properties"][0]
            smiles = (props.get("IsomericSMILES") or props.get("CanonicalSMILES")
                      or props.get("ConnectivitySMILES"))
            if not smiles:
                return {"similar_cids": [], "note": "no_smiles"}
        except (requests.RequestException, KeyError, IndexError, ValueError):
            return {"similar_cids": [], "note": "no_smiles"}
        try:
            r = requests.get(
                f"{self.base}/compound/fastsimilarity_2d/smiles/"
                f"{requests.utils.quote(smiles, safe='')}/cids/JSON",
                params={"Threshold": threshold, "MaxRecords": max_records},
                timeout=self.timeout)
            if r.status_code != 200:
                return {"similar_cids": [], "note": "no_similar"}
            return {"similar_cids": r.json()["IdentifierList"]["CID"]}
        except (requests.RequestException, KeyError, ValueError):
            return {"similar_cids": [], "note": "api_failure"}
