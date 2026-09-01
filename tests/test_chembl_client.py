"""ChEMBL API client tests without live network access."""
from ethno_evidence.data.chembl_client import ChEMBLClient


class _Response:
    def __init__(self, payload):
        self._payload = payload

    def raise_for_status(self):
        return None

    def json(self):
        return self._payload


class _Session:
    def __init__(self, payload):
        self.payload = payload
        self.calls = []

    def get(self, url, params, timeout):
        self.calls.append({"url": url, "params": params, "timeout": timeout})
        return _Response(self.payload)


def test_target_activities_applies_temporal_filters_and_keeps_provenance():
    payload = {
        "activities": [{
            "target_chembl_id": "CHEMBL203",
            "target_pref_name": "Cyclooxygenase-2",
            "standard_type": "IC50",
            "standard_value": "120",
            "standard_units": "nM",
            "pchembl_value": "6.92",
            "document_chembl_id": "CHEMBL112233",
            "document_year": 2019,
            "assay_chembl_id": "CHEMBL445566",
            "assay_type": "B",
        }]
    }
    client = ChEMBLClient()
    client._session = _Session(payload)

    out = client.target_activities(
        "CHEMBL25",
        organism="Homo sapiens",
        order_by="-pchembl_value",
        document_year_lte=2022,
        document_year_gte=2000,
    )

    params = client._session.calls[0]["params"]
    assert params["document_year__lte"] == 2022
    assert params["document_year__gte"] == 2000
    assert params["target_organism"] == "Homo sapiens"
    assert out[0]["pchembl_value"] == 6.92
    assert out[0]["document_chembl_id"] == "CHEMBL112233"
    assert out[0]["document_year"] == 2019
    assert out[0]["assay_chembl_id"] == "CHEMBL445566"
    assert out[0]["assay_type"] == "B"


def test_resolve_compound_targets_forwards_temporal_filters(monkeypatch):
    client = ChEMBLClient()
    monkeypatch.setattr(
        client,
        "search_molecules",
        lambda name, limit=1: [{"molecule_chembl_id": "CHEMBL1"}],
    )
    captured = {}

    def fake_target_activities(chembl_id, **kwargs):
        captured.update(kwargs)
        return [{"target_pref_name": "AKT1", "pchembl_value": 7.0}]

    monkeypatch.setattr(client, "target_activities", fake_target_activities)
    out = client.resolve_compound_targets(
        ["example"],
        sleep_s=0,
        document_year_lte=2022,
        document_year_gte=1990,
    )

    assert captured["document_year_lte"] == 2022
    assert captured["document_year_gte"] == 1990
    assert out["example"]["targets"] == ["AKT1"]
