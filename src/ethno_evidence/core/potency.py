"""数据源无关的效价换算与证据分桶。

论文口径（与 ChEMBL pchembl_value 阈值对应）：
- strong 档：效价 <= 10 uM   （等价 ChEMBL pchembl_value >= 5.0）
- weak 档：10 uM < 效价 <= 100 uM（等价 ChEMBL pchembl_value 4.0-5.0）

所有数据源（ChEMBL / TTD / 未来其他）通过换算函数把原始效价统一为
µM 线性值后，共用本模块分桶，保证跨源口径一致、结果可分桶对比。
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

# 证据分桶边界（µM），对应论文 "约10µM / 10–100µM" 口径
STRONG_UM = 10.0
WEAK_UM = 100.0

# 化合物-靶点活性类型（含被排除类型，用于审计记录）
ACTIVITY_TYPES = {"IC50", "Ki", "EC50"}
_ACTIVITY_TYPE_UPPER = {t.upper() for t in ACTIVITY_TYPES}


@dataclass
class PotencyRecord:
    """一条化合物-靶点实测活性记录，效价已统一为 µM。"""
    target: str              # 基因名（HUGO 习惯大写，如 RELA）
    potency_um: float        # 统一后的 µM 线性值
    activity_type: str       # IC50 / Ki / EC50
    raw: str                 # 原始字符串（如 "IC50 = 550 nM"），用于审计


def pchembl_to_um(pchembl_value: float) -> float:
    """ChEMBL pchembl_value (log10 M) -> µM。pchembl 5.0 => 10 µM。"""
    return 10.0 ** (6.0 - pchembl_value)


def nm_to_um(value_nm: float) -> float:
    """线性 nM -> µM。"""
    return value_nm / 1000.0


def parse_ttd_activity(raw: str) -> Optional[PotencyRecord]:
    """解析 TTD 的 Activity 列，如 'IC50 = 550 nM' / 'Ki < 10 nM'。

    保守规则：
    - 仅接受 nM 单位的 IC50 / Ki / EC50（µg/mL 等非摩尔单位无法换算，排除）；
    - '<' 前缀表示真实值更强，按边界值（最弱估计）参与分桶，安全；
    - '>' 前缀表示真实值更弱，按边界值会高估活性，排除；
    - 返回的 potency_um 统一为 µM。
    """
    if not raw or "=" not in raw:
        return None
    left, right = raw.split("=", 1)
    act_type = left.strip()
    if act_type.upper() not in _ACTIVITY_TYPE_UPPER:
        return None
    value_txt = right.strip()
    qualifier = ""
    if value_txt.startswith("<"):
        qualifier = "<"
        value_txt = value_txt[1:].strip()
    elif value_txt.startswith(">"):
        return None  # 下限不可靠，保守排除
    parts = value_txt.split()
    if len(parts) != 2 or parts[1].upper() != "NM":
        return None
    try:
        value_nm = float(parts[0])
    except ValueError:
        return None
    return PotencyRecord(
        target="",
        potency_um=nm_to_um(value_nm),
        activity_type=act_type,
        raw=f"{act_type} {qualifier}{value_txt}",
    )


def bucket_by_potency(records: List[PotencyRecord],
                      strong_um: float = STRONG_UM,
                      weak_um: float = WEAK_UM,
                      top_n: int = 30) -> Tuple[List[str], List[str]]:
    """按效价把靶点分为 strong / weak 两档（论文口径）。

    - 每个靶点取多条记录中的最优（最小 µM）值，避免单条高噪记录拉低档位；
    - strong 档按效价升序截断到 top_n，与 ChEMBL 路径行为一致；
    - weak 档不截断，仅供敏感性分析。
    返回 (strong_targets, weak_targets)，两档分开、互不混入。
    """
    best: Dict[str, float] = {}
    for rec in records:
        if rec.potency_um is None or not rec.target:
            continue
        if rec.target not in best or rec.potency_um < best[rec.target]:
            best[rec.target] = rec.potency_um
    strong = sorted((t for t, um in best.items() if um <= strong_um),
                    key=lambda t: best[t])[:top_n]
    weak = sorted((t for t, um in best.items() if strong_um < um <= weak_um),
                  key=lambda t: best[t])
    return strong, weak
