"""Validate the ingested Torabi & Kruglyak 2012 tables against the paper and themselves.

    python3 -I -m pytest -q tests/test_torabi2012.py      (after ingest/torabi2012.py)
"""
import json
import os
import re

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "sources", "torabi2012")
R64 = {"chrI": 230218, "chrII": 813184, "chrIII": 316620, "chrIV": 1531933, "chrV": 576874, "chrVI": 270161, "chrVII": 1090940,
       "chrVIII": 562643, "chrIX": 439888, "chrX": 745751, "chrXI": 666816, "chrXII": 1078177, "chrXIII": 924431, "chrXIV": 784333,
       "chrXV": 1091291, "chrXVI": 948066}


def rows(tab):
    keys = list(tab)
    return [dict(zip(keys, v)) for v in zip(*(tab[k] for k in keys))]


def tsv(name):
    with open(os.path.join(SRC, "extracted", name)) as f:
        lines = [l.rstrip("\n").split("\t") for l in f if not l.startswith("#")]
    return lines[0], lines[1:]


@pytest.fixture(scope="module")
def D():
    with open(os.path.join(ROOT, "docs", "data.json"), encoding="utf-8") as f:
        return json.load(f)


@pytest.fixture(scope="module")
def ds(D):
    return D["datasets"]["id"].index("torabi2012")


@pytest.fixture(scope="module")
def Q(D, ds):
    return [r for r in rows(D["qtl"]) if r["ds"] == ds]


def page_text():
    t = open(os.path.join(SRC, "article.html"), encoding="utf-8", errors="replace").read()
    import html
    return re.sub(r"\s+", " ", html.unescape(re.sub(r"<[^>]+>", " ", t)))


def test_quotes_verbatim():
    head, fs = tsv("text_facts.tsv")
    text = page_text()
    assert len(fs) == 15
    for r in fs:
        assert dict(zip(head, r))["quote"] in text


def test_strains_and_two_crosses(D, ds):
    X = [r for r in rows(D["crosses"]) if r["ds"] == ds]
    assert {x["code"] for x in X} == {"WT", "sup35"}
    P = D["parents"]
    assert all({P["name"][x["p1"]], P["name"][x["p2"]]} == {"BY", "RM"} for x in X)
    assert all(x["n"] is None for x in X)


def test_nine_xqtl_conditions(D, ds):
    head, s1 = tsv("table_s1.tsv")
    assert len(s1) == 32                                   # conditions tested in Table S1
    mapped = [r for r in s1 if r[head.index("X-QTL Concentration")] not in ("", "NA")]
    assert len(mapped) == 9                                # "the nine conditions"
    assert sum(1 for t in D["traits"]["ds"] if t == ds) == 9


def test_chrXIII_locus_rows(D, Q):
    assert len(Q) == 2
    T, P = D["traits"], D["parents"]
    assert {T["raw"][q["t"]].split(";")[0] for q in Q} == {"Diamide", "H2O2"}
    for q in Q:
        assert q["c"] == "chrXIII" and q["a"] is None and q["e"] is None and q["p"] is None
        assert q["m"] == "x" and P["name"][q["up"]] == "RM"
        assert D["crosses"]["code"][q["x"]] == "sup35"


def test_sky1_tier_x_on_chrXIII(D, ds, Q):
    G = D["genes"]
    for q in Q:
        assert G["n"][q["g"]] == "SKY1" and q["gt"] == "x" and G["c"][q["g"]] == "chrXIII"
    C = [r for r in rows(D["causal"]) if r["ds"] == ds]
    assert len(C) == 2 and all(c["tier"] == "x" and c["name"] == "SKY1" for c in C)


def test_direction_matches_table_s2():
    head, s2 = tsv("table_s2.tsv")
    est = {r[0]: float(r[1]) for r in s2}
    assert est["SKY1 (BY)"] < 0          # the BY allele of SKY1 lowers growth on H2O2, i.e. RM raises it
    assert est["BG (BY)"] < 0            # "RM growing better ... than BY"


@pytest.mark.xfail(strict=True, reason=(
    "Text: 18 readthrough-dependent locus-condition pairs, ten distinct loci ('X-QTL reveals readthrough-dependent loci'). "
    "Only the chromosome XIII locus (diamide, H2O2) is stated in the text; the rest are drawn in Figure 2 and Figure S2 "
    "(raw/journal.pgen.1002546.s002.tif) and not tabulated, so the atlas holds 2 rows."))
def test_18_readthrough_dependent_loci(Q):
    assert len(Q) == 18


@pytest.mark.skip(reason="'one to 20' loci per condition per cross are drawn in Figure S2 only; not checkable from tables")
def test_loci_per_condition():
    pass


def test_provenance_and_verbatim(D, ds):
    codes = set(D["sources"]["code"])
    sc = D["datasets"]["sc"][ds]
    kept = {}
    for ch in R64:
        with open(os.path.join(ROOT, "docs", "v", f"{sc}_{ch}.json"), encoding="utf-8") as f:
            v = json.load(f)
        for code, rs in v["rows"].items():
            kept.update(rs)
    for t in ("qtl", "causal"):
        for r in rows(D[t]):
            if r["ds"] != ds:
                continue
            code, row = r["pv"].split("#")
            assert code in codes and code.startswith("t12:") and int(row) >= 2
            assert r["id"] in kept
    with open(os.path.join(ROOT, "docs", "v", f"{sc}_extra.json"), encoding="utf-8") as f:
        ex = json.load(f)
    assert len(ex["t12:S1"]["rows"]) == 32 and len(ex["t12:S2"]["rows"]) == 8 and len(ex["t12:S3"]["rows"]) == 8
