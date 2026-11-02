"""Validate the Brem and Kruglyak 2005 (PNAS) ingest against the paper's own text and tables.

    python3 -I -m pytest -q tests/test_brem2005pnas.py      (after ./build.sh)

The article reports summaries only, so the checks are internal consistency between the numbers in
the text, the abstract and Tables 1 and 2, plus provenance. Disagreements are strict xfails naming
the source and column; the data are not changed to match.
"""
import csv
import html
import json
import os
import re

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "sources", "brem2005pnas")


@pytest.fixture(scope="module")
def D():
    with open(os.path.join(ROOT, "docs", "data.json"), encoding="utf-8") as f:
        return json.load(f)


@pytest.fixture(scope="module")
def ds(D):
    return D["datasets"]["id"].index("brem2005pnas")


@pytest.fixture(scope="module")
def facts():
    with open(os.path.join(SRC, "extracted", "text_facts.tsv"), encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


@pytest.fixture(scope="module")
def tables():
    out = {}
    with open(os.path.join(SRC, "extracted", "tables.tsv"), encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            out.setdefault(r["table"], []).append([r[f"c{i}"] for i in range(1, 6)])
    return out


def num(facts, start):
    v = next(d["value"] for d in facts if d["value"].startswith(start))
    return int(re.match(r"[\d,]+", v).group(0).replace(",", ""))


def test_every_quote_is_verbatim_in_the_article(facts):
    t = open(os.path.join(SRC, "article.html"), encoding="utf-8").read()
    body = re.sub(r"<(script|style)[^>]*>.*?</\1>", " ", t, flags=re.S)
    txt = re.sub(r"\s+", " ", html.unescape(re.sub(r"<[^>]+>", " ", body)))
    assert len(facts) == 20
    for d in facts:
        assert d["quote"] in txt, d["value"]


def test_design(D, ds):
    X = [i for i, d in enumerate(D["crosses"]["ds"]) if d == ds]
    assert len(X) == 1 and D["crosses"]["n"][X[0]] == 112
    P = D["parents"]
    assert sorted(P["name"][i] for i, d in enumerate(P["ds"]) if d == ds) == ["BY", "RM"]
    assert D["datasets"]["fdr"][ds] == 0.05


def test_no_qtl_rows_are_invented(D, ds):
    assert D["datasets"]["nq"][ds] == 0
    assert ds not in D["qtl"]["ds"]
    assert not D.get("lazy", {}).get("ds") or ds not in D["lazy"]["ds"]


def test_heritable_counts_match_percentages(facts):
    h = num(facts, "3,546 highly heritable")
    for start, pct in (("2,091", 59), ("406", 11), ("2,093", 59), ("583", 16)):
        assert round(100 * num(facts, start) / h) == pct, start


@pytest.mark.xfail(strict=True, reason=(
    "Abstract: no QTL for '40%' of highly heritable transcripts. Results: 2,091 of 3,546 (59%) have a QTL "
    "(extracted/text_facts.tsv, rows '3,546 ...' and '2,091 ...'), so 1,455 (41.0%) have none."))
def test_no_qtl_fraction_is_40_percent(facts):
    h, w = num(facts, "3,546 highly heritable"), num(facts, "2,091")
    assert round(100 * (h - w) / h) == 40


def test_table1_supports_text(tables):
    t1 = {r[0]: r for r in tables["Table 1"][1:]}
    data = {k: float(v[1]) for k, v in t1.items() if k != "Median"}
    assert abs(sum(data.values()) - 1) < 0.011
    assert float(t1["Median"][1]) == 0.27                     # median variance explained 27%
    assert round(data["50–60"] + data[">60"], 2) == 0.23     # strongest QTL > 50% for 23%
    assert data["0–10"] + data["10–20"] + data["20–30"] > 0.5  # < 30% for more than half


def test_table2_supports_abstract(tables):
    t2 = {int(r[0]): [float(x) for x in r[1:5]] for r in tables["Table 2"][1:]}
    assert t2[1][0] == 0.03                                   # 3% single locus
    assert round(t2[1][3] + t2[2][3], 2) == 0.18              # at most 17-18% one or two loci
    assert abs(t2[5][1] - 0.5) < 0.02                         # half need more than five loci
    assert abs(t2[8][1] - 0.30) < 0.005                       # almost a third need more than eight


def test_provenance_and_extra(D, ds):
    codes = {c for c, d in zip(D["sources"]["code"], D["sources"]["ds"]) if d == ds}
    assert codes == {"b05p:TF", "b05p:T", "b05p:JATS"}
    for pv, d in zip(D["crosses"]["pv"], D["crosses"]["ds"]):
        if d == ds:
            assert pv.split("#")[0] in codes
    with open(os.path.join(ROOT, "docs", "v", "b05p_extra.json"), encoding="utf-8") as f:
        ex = json.load(f)
    assert len(ex["b05p:TF"]["rows"]) == 20 and len(ex["b05p:T"]["rows"]) == 20
