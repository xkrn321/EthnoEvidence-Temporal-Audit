"""TTD (Therapeutic Target Database) 本地 TSV 客户端（只读、无需网络/密钥）。

数据来源：TTD 官网免费下载文件（https://ttd.idrblab.cn/download，Version 10.1.01）：
- P1-03-TTD_crossmatching.txt  药物 ID <-> 名称 / PubChem CID 等
- P1-01-TTD_target_download.txt 靶点 ID <-> 基因名 / UniProt 等
- P1-09-Target_compound_activity.txt 化合物-靶点活性（IC50/Ki/EC50, nM）

三个文件均为长格式（ID\\t字段名\\t值）。解析后提供与 ChEMBL 路径同构的
resolve_compound_targets()：成分名 -> 靶点基因名列表（按 strong/weak 两档
分桶，共用 ethno_evidence.core.potency 的统一口径）。
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, List, Optional

from ..core.potency import PotencyRecord, bucket_by_potency, parse_ttd_activity

DEFAULT_DATA_DIR = Path(__file__).resolve().parents[3] / "data" / "ttd"

# 名称规范化：
# 1) 去掉开头的数字+连字符前缀（如 "6-gingerol" -> "gingerol"）；
# 2) 去掉希腊字母前缀（如 "beta-sitosterol" -> "sitosterol"）。
# 压缩匹配：去空格/下划线/连字符/括号后比较（如 "oleanolic_acid" <-> "oleanolic acid"）。
# 均仅作为精确匹配失败后的第二、三层匹配，防止误配。
_NUM_PREFIX_RE = re.compile(r"^[0-9]+[-,]?")
_GREEK_PREFIX_RE = re.compile(r"^(?:beta|alpha|gamma|delta)-?")


def _normalize_name(name: str) -> str:
    n = name.strip().lower()
    m = _NUM_PREFIX_RE.match(n)
    if m:
        n = n[m.end():].lstrip("- ")
    m = _GREEK_PREFIX_RE.match(n)
    if m:
        n = n[m.end():].lstrip("- ")
    return n


def _compact(name: str) -> str:
    return re.sub(r"[\s_\-\(\)]", "", name.strip().lower())


class TTDClient:
    """加载 TTD 免费下载 TSV 并解析化合物-靶点活性。"""

    def __init__(self, data_dir=DEFAULT_DATA_DIR):
        self.data_dir = Path(data_dir)
        # 索引 1：化合物 ID -> 名称/PubChem CID（来自 P1-03）
        self.comp_name: Dict[str, str] = {}          # ttd_drug_id -> 规范名
        self.comp_cid: Dict[str, str] = {}           # ttd_drug_id -> pubchem cid
        self.name_to_id: Dict[str, str] = {}         # 小写名称 -> ttd_drug_id
        self.name_to_id_compact: Dict[str, str] = {}  # 压缩名称 -> ttd_drug_id
        # 索引 2：靶点 ID -> 基因名（来自 P1-01）
        self.target_gene: Dict[str, str] = {}
        # 索引 3：活性记录（来自 P1-09）
        self.activities_by_id: Dict[str, List[str]] = {}    # ttd_drug_id -> [raw activity]
        self.activities_by_cid: Dict[str, List[str]] = {}   # pubchem cid -> [raw activity]
        self._loaded = False

    # ------------------------------------------------------------------ load
    def load(self) -> None:
        if self._loaded:
            return
        self._load_crossmatching()
        self._load_targets()
        self._load_activities()
        self._loaded = True

    def _load_crossmatching(self) -> None:
        path = self.data_dir / "P1-03-TTD_crossmatching.txt"
        if not path.exists():
            raise FileNotFoundError(f"TTD crossmatching 文件缺失: {path}")
        with path.open(encoding="utf-8", errors="ignore") as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) != 3:
                    continue
                drug_id, field, value = parts
                if not drug_id.startswith("D") or not value.strip():
                    continue
                if field == "DRUGNAME":
                    self.comp_name.setdefault(drug_id, value.strip())
                    key = value.strip().lower()
                    self.name_to_id.setdefault(key, drug_id)
                    self.name_to_id_compact.setdefault(_compact(key), drug_id)
                elif field == "PUBCHCID":
                    self.comp_cid.setdefault(drug_id, value.strip())

    def _load_targets(self) -> None:
        path = self.data_dir / "P1-01-TTD_target_download.txt"
        if not path.exists():
            raise FileNotFoundError(f"TTD target 文件缺失: {path}")
        with path.open(encoding="utf-8", errors="ignore") as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 3 and parts[0].startswith("T") and parts[1] == "GENENAME":
                    if parts[2].strip():
                        self.target_gene[parts[0]] = parts[2].strip()

    def _load_activities(self) -> None:
        path = self.data_dir / "P1-09-Target_compound_activity.txt"
        if not path.exists():
            raise FileNotFoundError(f"TTD activity 文件缺失: {path}")
        with path.open(encoding="utf-8", errors="ignore") as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) != 4:
                    continue
                target_id, comp_id, cid, activity = parts
                if not target_id.startswith("T") or "=" not in activity:
                    continue
                self.activities_by_id.setdefault(comp_id, []).append((target_id, activity))
                if cid.strip():
                    self.activities_by_cid.setdefault(cid.strip(), []).append((target_id, activity))

    # ------------------------------------------------------------- resolution
    def resolve_compound_targets(self, compound_names: List[str],
                                 top_n: int = 30) -> dict:
        """成分名 -> 靶点基因名（strong/weak 两档）。

        返回结构与 ChEMBL 路径同构：
        {name: {"ttd_id": str|None, "records": [PotencyRecord...],
                "strong": [...], "weak": [...], "error"?: str}}
        名称未命中/无活性记录时记录为 None/空，绝不虚构。
        """
        self.load()
        result: Dict[str, dict] = {}
        for name in compound_names:
            entry = {"ttd_id": None, "records": [], "strong": [], "weak": []}
            key = name.strip().lower()
            drug_id = self.name_to_id.get(key)
            if not drug_id:
                drug_id = self.name_to_id_compact.get(_compact(key))
            if not drug_id:  # 去数字/希腊字母前缀后的规范化名称再尝试
                norm = _normalize_name(name)
                drug_id = self.name_to_id.get(norm) or self.name_to_id_compact.get(_compact(norm))
            raw_acts: List = []
            if drug_id:
                entry["ttd_id"] = drug_id
                raw_acts = self.activities_by_id.get(drug_id, [])
            if not raw_acts and drug_id:
                cid = self.comp_cid.get(drug_id)
                if cid:
                    raw_acts = self.activities_by_cid.get(cid, [])
            if not drug_id:
                entry["error"] = "name_not_found"
            elif not raw_acts:
                entry["error"] = "no_activity"
            records: List[PotencyRecord] = []
            for target_id, activity in raw_acts:
                parsed = parse_ttd_activity(activity)
                if parsed is None:
                    continue
                gene = self.target_gene.get(target_id)
                if not gene:
                    continue
                parsed.target = gene
                records.append(parsed)
            entry["records"] = records
            entry["strong"], entry["weak"] = bucket_by_potency(records, top_n=top_n)
            result[name] = entry
        return result
