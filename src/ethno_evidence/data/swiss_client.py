"""SwissTargetPrediction (STP) 网页服务客户端（predicted 等级证据）。

STP 无官方 API；本客户端模拟浏览器表单提交流程：
1) GET 首页获取 session cookie；
2) POST /predict.php（smiles + organism + ioi=2）提交预测任务；
3) 从等待页提取结果 URL（/result.php?job=...&organism=...）；
4) 轮询结果页直到计算完成（出现 #resultTable，STP 提示最多约 1 分钟）；
5) 解析靶点表（Target / Common name / Uniprot / ChEMBL / Target Class /
   Probability / Known actives 3D/2D）。

证据等级：STP 基于 2D/3D 结构相似性（SEA 原理）预测靶点，属 predicted 层，
严格与 measured（ChEMBL/TTD 实测效价）分档，绝不混入。
"""
from __future__ import annotations

import re
import time
from concurrent.futures import ThreadPoolExecutor
from typing import Dict, List, Optional

import requests
from bs4 import BeautifulSoup

BASE = "https://www.swisstargetprediction.ch"
_UA = ("Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 "
       "(KHTML, like Gecko) Chrome/125.0 Safari/537.36")


def parse_targets(html: str) -> List[dict]:
    """解析 STP 结果页中的 #resultTable 为靶点列表（无则返回空）。"""
    soup = BeautifulSoup(html, "html.parser")
    table = soup.find("table", id="resultTable")
    if not table:
        return []
    out = []
    for tr in table.find_all("tr"):
        tds = tr.find_all("td")
        if len(tds) < 7:
            continue
        out.append({
            "target": tds[0].get_text(strip=True),
            "common_name": tds[1].get_text(strip=True),
            "uniprot": tds[2].get_text(strip=True),
            "chembl_id": tds[3].get_text(strip=True),
            "target_class": tds[4].get_text(strip=True),
            "probability": _as_float(tds[5].get_text(strip=True)),
            "known_actives": tds[6].get_text(strip=True),
        })
    return out


def _as_float(value) -> Optional[float]:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


class SwissTargetPredictionClient:
    def __init__(self, timeout_s: int = 30, poll_interval_s: float = 2.0,
                 max_wait_s: float = 150.0, offline: bool = False):
        self.timeout = timeout_s
        self.poll_interval = poll_interval_s
        self.max_wait = max_wait_s
        self.offline = offline

    def predict(self, smiles: str, organism: str = "Homo_sapiens") -> dict:
        """对单个 SMILES 提交 STP 预测并轮询取回结果。

        结果按 SMILES 磁盘缓存（configs/cache/stp_*），重复分析秒回、
        不重击 STP 服务器。返回 {"job": int|None, "targets": [...],
        "error": str|None}；失败时 targets 为空并给出原因，绝不虚构。
        """
        if self.offline:
            return {"job": None, "targets": [], "error": "offline"}
        from .cache import cache_get, cache_set
        key = f"{smiles}|{organism}"
        cached = cache_get("stp", key)
        if cached is not None:
            return cached
        result = self._predict_remote(smiles, organism=organism)
        cache_set("stp", key, result)
        return result

    def _predict_remote(self, smiles: str, organism: str = "Homo_sapiens") -> dict:
        """实际在线提交（无缓存路径）。"""
        session = requests.Session()
        session.headers.update({"User-Agent": _UA})
        # 1) 取 session cookie
        try:
            session.get(BASE + "/", timeout=self.timeout)
            # 2) 提交任务（ioi=2 是页面 formSubmit() 写入的必需字段）
            r = session.post(
                BASE + "/predict.php",
                data={"smiles": smiles, "organism": organism, "ioi": "2"},
                headers={"Referer": BASE + "/"},
                timeout=self.timeout)
            r.raise_for_status()
        except requests.RequestException:
            return {"job": None, "targets": [], "error": "submit_failed"}
        m = re.search(r"(/result\.php\?job=\d+&organism=[A-Za-z_]+)", r.text)
        if not m:
            return {"job": None, "targets": [], "error": "no_job_url"}
        job_url = BASE + m.group(1)
        job_id = int(re.search(r"job=(\d+)", job_url).group(1))
        # 3) 轮询结果页
        deadline = time.time() + self.max_wait
        last_html = ""
        while time.time() < deadline:
            time.sleep(self.poll_interval)
            try:
                rr = session.get(job_url, timeout=self.timeout)
                if rr.status_code != 200:
                    continue
            except requests.RequestException:
                continue
            last_html = rr.text
            if "resultTable" in rr.text:
                return {"job": job_id, "targets": parse_targets(rr.text), "error": None}
            # 结果已过期/被清理
            if "not available" in rr.text or "has expired" in rr.text.lower():
                return {"job": job_id, "targets": [], "error": "job_expired"}
        return {"job": job_id, "targets": [], "error": "timeout"}

    def predict_many(self, name_smiles: Dict[str, str],
                     organism: str = "Homo_sapiens",
                     per_request_delay_s: float = 0.8,
                     workers: int = 4) -> Dict[str, dict]:
        """对多个成分（name -> SMILES）并发预测；保持字典序输出。

        命中缓存的成分直接返回；未命中的并发提交（STP 单任务上限约 1 分钟）。
        返回 {name: {"job", "targets", "error"}}。
        """
        # 先剥出缓存命中项，避免重复提交
        from .cache import cache_get
        out: Dict[str, dict] = {}
        pending: Dict[str, str] = {}
        for name, smiles in name_smiles.items():
            cached = cache_get("stp", f"{smiles}|{organism}")
            if cached is not None:
                out[name] = cached
            else:
                pending[name] = smiles

        def _run(item):
            n, s = item
            return n, self.predict(s, organism=organism)

        with ThreadPoolExecutor(max_workers=workers) as ex:
            for n, res in ex.map(_run, list(pending.items())):
                out[n] = res
                if per_request_delay_s > 0:
                    time.sleep(per_request_delay_s)
        # 保持输入顺序
        return {n: out[n] for n in name_smiles if n in out}
