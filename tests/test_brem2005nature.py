"""Validate the ingested Brem et al. 2005 (Nature) rows against the paper's own statements.

    python3 -I -m pytest -q tests/test_brem2005nature.py      (after the ingest)
"""
import csv
import html
import json
import math
import os
import re

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "sources", "brem2005nature")


def rows(tab):
    keys = list(tab)
    return [dict(zip(keys, v)) for v in zip(*(tab[k] for k in keys))]


@pytest.fixture(scope="module")
def D():
    with open(os.path.join(ROOT, "docs", "data.json"), encoding="utf-8") as f:
        return json.load(f)


@pytest.fixture(scope="module")
def ds(D):
    return D["datasets"]["id"].index("brem2005nature")


@pytest.fixture(scope="module")
def facts():
    with open(os.path.join(SRC, "extracted", "text_facts.tsv")) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def page(path):
    t = open(path, encoding="utf-8", errors="replace").read()
    t = re.sub(r"<(script|style)[^>]*>.*?</\1>", " ", t, flags=re.S)
    return re.sub(r"\s+", " ", html.unescape(re.sub(r"<[^>]+>", " ", t)))


def num(facts, kind, pattern):
    for f in facts:
        if f["kind"] == kind:
            m = re.search(pattern, f["quote"])
            if m:
                return [float(g) for g in m.groups()]
    raise AssertionError(pattern)


def test_quotes_verbatim(facts):
    pages = {"raw/PMC1409747.html": page(os.path.join(SRC, "raw", "PMC1409747.html")), "article.html": page(os.path.join(SRC, "article.html"))}
    assert len(facts) == 26
    for f in facts:
        assert f["quote"] in pages[f["file"]], f["quote"][:60]


def test_design(D, ds):
    X = [r for r in rows(D["crosses"]) if r["ds"] == ds]
    assert len(X) == 1 and X[0]["n"] == 112
    P = [r for r in rows(D["parents"]) if r["ds"] == ds]
    assert sorted(p["name"] for p in P) == ["BY", "RM"]
    assert D["datasets"]["fdr"][ds] == 0.05


def test_counts_consistent(facts):
    n225, = num(facts, "method", r"(\d+) expression traits were called significant")
    a, b, c = num(facts, "secondary", r"only (\d+) \((\d+)%\) of the (\d+) were")
    assert n225 == c == 225 and a == 74 and round(100 * a / c) == b == 33
    assert "67% of secondary loci" in next(f["quote"] for f in facts if f["kind"] == "secondary" and "74" in f["quote"])
    assert 0.05 * 225 < 12   # "fewer than 12 of these are expected to have either locus as a false positive"


WORDS = {"two": 2, "four": 4, "seven": 7}


def wnum(v):
    return float(v) if v[0].isdigit() else float(WORDS[v])


def test_mat_gpa1_pair_counts(facts):
    q = {f["quote"][:30]: f["quote"] for f in facts if f["kind"] == "pair"}
    big, extra = [wnum(v) for v in re.search(r"was (\d+), with (\w+) additional", " ".join(q.values())).groups()]
    assert big == 14 and big + extra - 1 == 17   # 18 transcripts, minus the silenced HMLα1 copy
    m = re.search(r"(\d+) were previously shown to be regulated by mating type, of which (\w+) are also regulated by the "
                  r"pheromone response pathway; (\w+) have known function that does not involve mating, and (\w+) have no known function",
                  " ".join(q.values()))
    known, pher, other, unknown = [wnum(v) for v in m.groups()]
    assert (known, pher, other, unknown) == (11, 7, 2, 4) and known + other + unknown == 17
    e1, e2, e3, e4 = num(facts, "engineered", r"For (\d+) of (\d+) transcripts previously known to be regulated by mating type, and (\d+) of (\d+)")
    assert (e1, e2, e3, e4) == (7, 11, 3, 4) and e2 == known and e4 == unknown


def test_bins_match_r64(D, facts):
    nb, = num(facts, "bins", r"divided into (\d+) bins of 20 kb")
    assert nb == 611 == sum(math.ceil(v / 20000) for v in D["chroms"].values())


def test_qtl_rows(D, ds):
    Q = [r for r in rows(D["qtl"]) if r["ds"] == ds]
    assert len(Q) == 2
    assert [q["c"] for q in Q] == ["chrIII", "chrVIII"]
    for q in Q:
        assert q["a"] is None and q["e"] is None and q["p"] is None
        assert q["cc"] == 17 and q["gt"] == "x" and q["m"] == "w"
        assert D["genes"]["c"][q["g"]] == q["c"]


def test_confirmed_genes(D, ds):
    C = [r for r in rows(D["causal"]) if r["ds"] == ds]
    assert sorted(c["orf"] for c in C) == ["YCR040W", "YHR005C"]
    assert all(c["tier"] == "x" for c in C)
    g = {c["orf"]: c for c in C}
    assert g["YHR005C"]["c"] == "chrVIII" and g["YCR040W"]["c"] == "chrIII"


@pytest.mark.xfail(strict=True, reason=(
    "The text names 'YOR090C/PTC1' (extracted/text_facts.tsv, kind engineered). In SGD R64-5-1 (sources/sgd/"
    "genes_R64-5-1.tsv) YOR090C is PTC5 and PTC1 is YDL006W; the paper's systematic and standard names disagree."))
def test_yor090c_is_ptc1(D):
    G = D["genes"]
    assert G["n"][G["id"].index("YOR090C")] == "PTC1"


def test_provenance(D, ds):
    codes = {c for c, d in zip(D["sources"]["code"], D["sources"]["ds"]) if d == ds}
    for t in ("crosses", "qtl", "causal"):
        for r in rows(D[t]):
            if r.get("ds", ds) != ds or not r.get("pv"):
                continue
            code, n = r["pv"].split("#")
            assert code in codes and int(n) >= 2


def test_verbatim_rows_kept(D, ds):
    ids = [r["id"] for r in rows(D["qtl"]) if r["ds"] == ds] + [r["id"] for r in rows(D["causal"]) if r["ds"] == ds]
    seen = set()
    for ch in D["chroms"]:
        with open(os.path.join(ROOT, "docs", "v", f"b05n_{ch}.json"), encoding="utf-8") as f:
            seen |= set(json.load(f)["rows"]["b05n:TF"])
    assert set(ids) == seen
    with open(os.path.join(ROOT, "docs", "v", "b05n_extra.json"), encoding="utf-8") as f:
        assert len(json.load(f)["b05n:TF"]["rows"]) == 26
