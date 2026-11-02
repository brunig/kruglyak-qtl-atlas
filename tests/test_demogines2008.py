"""Validate the ingested Demogines et al. 2008 tables against the article.

    python3 -I -m pytest -q tests/test_demogines2008.py      (after ./build.sh)

Numbers come from the article text (quoted in sources/demogines2008/extracted/text_facts.tsv).
Where the article disagrees with itself, the test is a strict xfail naming where.
"""
import html
import json
import os
import re

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "sources", "demogines2008")
R64 = {"chrI": 230218, "chrII": 813184, "chrIII": 316620, "chrIV": 1531933, "chrV": 576874, "chrVI": 270161,
       "chrVII": 1090940, "chrVIII": 562643, "chrIX": 439888, "chrX": 745751, "chrXI": 666816, "chrXII": 1078177,
       "chrXIII": 924431, "chrXIV": 784333, "chrXV": 1091291, "chrXVI": 948066}


def rows(tab):
    keys = list(tab)
    return [dict(zip(keys, v)) for v in zip(*(tab[k] for k in keys))]


@pytest.fixture(scope="module")
def D():
    with open(os.path.join(ROOT, "docs", "data.json"), encoding="utf-8") as f:
        return json.load(f)


@pytest.fixture(scope="module")
def ds(D):
    return D["datasets"]["id"].index("demogines2008")


@pytest.fixture(scope="module")
def Q(D, ds):
    return [r for r in rows(D["qtl"]) if r["ds"] == ds]


@pytest.fixture(scope="module")
def genes(D):
    g = D["genes"]
    return {g["n"][i] or g["id"][i]: (g["c"][i], g["a"][i], g["e"][i]) for i in range(len(g["id"]))}


def facts():
    with open(os.path.join(SRC, "extracted", "text_facts.tsv"), encoding="utf-8") as f:
        r = [l.rstrip("\n").split("\t") for l in f]
    return [dict(zip(r[0], x)) for x in r[1:]]


def article_flat():
    t = open(os.path.join(SRC, "article.html"), encoding="utf-8").read()
    t = re.sub(r"<(script|style)[^>]*>.*?</\1>", " ", t, flags=re.S)
    t = re.sub(r"</?(p|div|h[1-6]|li|tr|table|br|section)[^>]*>", "\n", t)
    t = html.unescape(re.sub(r"<[^>]+>", "", t))
    return re.sub(r"\s+", " ", t)


def test_quotes_are_verbatim():
    flat = article_flat()
    F = facts()
    assert len(F) == 31
    for f in F:
        assert flat.count(f["quote"]) == 1, f["quote"][:60]


def test_citation(D, ds):
    assert D["datasets"]["doi"][ds] == "10.1371/journal.pgen.1000123"
    page = open(os.path.join(SRC, "article.html"), encoding="utf-8").read()
    authors = re.findall(r'name="citation_author" content="([^"]+)"', page)
    assert authors == ["Ann Demogines", "Erin Smith", "Leonid Kruglyak", "Eric Alani"]


def test_design(D, ds):
    P = [r for r in rows(D["parents"]) if r["ds"] == ds]
    X = [r for r in rows(D["crosses"]) if r["ds"] == ds]
    assert [p["name"] for p in P] == ["BY", "RM"]
    assert {x["code"]: x["n"] for x in X} == {"BYxRM": 121, "BK3": 80}   # 121 genotyped segregants; 40 + 40 pooled


def test_traits(D, ds):
    T = [r for r in rows(D["traits"]) if r["ds"] == ds]
    assert sorted(t["ck"] for t in T) == ["4nqo 0.25ug/mL", "bleomycin 10ug/mL", "caffeine 10mM", "mms 0.03%"]
    assert all(t["cat"] == "drug/chemical" for t in T)


def test_six_qtl(Q):
    assert len(Q) == 6
    assert sum(q["m"] == "w" for q in Q) == 4 and sum(q["m"] == "x" for q in Q) == 2


def test_linkage_lods_above_thresholds(D, Q):
    tn = {i: D["traits"]["ck"][i] for i in range(len(D["traits"]["id"]))}
    thr = {"mms 0.03%": 3.9, "4nqo 0.25ug/mL": 3.8, "bleomycin 10ug/mL": 3.8, "caffeine 10mM": 3.7}
    lod = {tn[q["t"]]: q["lod"] for q in Q if q["m"] == "w"}
    assert lod == {"mms 0.03%": 16.6, "4nqo 0.25ug/mL": 10.1, "bleomycin 10ug/mL": 5.2, "caffeine 10mM": 5.0}
    assert all(lod[k] > thr[k] for k in lod)


def test_mms_and_4nqo_share_the_peak_marker(Q):
    W = [q for q in Q if q["lv"] == "YLR034C_1989"]
    assert len(W) == 2 and W[0]["c"] == W[1]["c"] == "chrXII" and W[0]["p"] == W[1]["p"]


def test_chrII_loci_do_not_overlap(Q):
    a, b = sorted((q for q in Q if q["c"] == "chrII"), key=lambda q: q["a"])
    assert a["e"] < b["a"]   # "their support intervals do not overlap"


def test_peak_reading_matches_stated_distances(Q, genes):
    p = next(q["p"] for q in Q if q["lv"] == "YLR034C_1989")
    _, ra, re_ = genes["RAD5"]
    _, aa, ae = genes["AAT2"]
    d_rad5 = min(abs(p - ra), abs(p - re_))
    d_aat2 = min(abs(p - aa), abs(p - ae))
    assert 1000 <= d_rad5 <= 3000      # "RAD5, located ~2 kb away from the peak marker"
    assert 11000 <= d_aat2 <= 14000    # "approximately 12.5 kb from the peak marker"


def test_bulk_segregant_regions(Q, genes):
    B = sorted((q for q in Q if q["m"] == "x"), key=lambda q: q["c"])
    assert [(q["c"], q["a"], q["p"], q["e"]) for q in B] == [("chrVII", 433546, 463019, 551683), ("chrXIV", 437935, 467209, 473640)]
    assert round((B[0]["e"] - B[0]["a"]) / 1000) == 118 and round((B[1]["e"] - B[1]["a"]) / 1000) == 36
    for q, g in zip(B, ("KAP122", "MKT1")):   # "the peak includes KAP122"; "the peak lies near MKT1" (R64-5-1)
        c, a, e = genes[g]
        assert c == q["c"] and a <= q["p"] <= e


def test_coordinates_inside_r64(Q):
    for q in Q:
        assert 1 <= q["a"] <= q["p"] <= q["e"] <= R64[q["c"]], q["id"]


def test_candidates_and_direction(D, ds, Q):
    G = D["genes"]
    tn = {i: D["traits"]["ck"][i] for i in range(len(D["traits"]["id"]))}
    P = D["parents"]
    got = {(G["n"][q["g"]], tn[q["t"]], q["m"], P["name"][q["up"]] if q["up"] >= 0 else None, q["gt"]) for q in Q if q["g"] >= 0}
    assert got == {("RAD5", "mms 0.03%", "w", "BY", "x"), ("RAD5", "4nqo 0.25ug/mL", "w", "BY", "x"),
                   ("MKT1", "4nqo 0.25ug/mL", "x", "RM", "x")}
    C = [r for r in rows(D["causal"]) if r["ds"] == ds]
    assert len(C) == 3 and {c["tier"] for c in C} == {"x"}


@pytest.mark.xfail(strict=True, reason=(
    "The Figure 2 legend says the MMS- and 4-NQO-sensitive phenotypes are 'linked to a locus on chromosome 11'; "
    "the Results text says chromosome 12 near SMF3 (marker YLR034C_1989), and SMF3 and RAD5 are on chrXII. "
    "The atlas follows the text."))
def test_figure2_legend_chromosome(Q):
    assert next(q["c"] for q in Q if q["lv"] == "YLR034C_1989") == "chrXI"


@pytest.mark.skip(reason="'Twenty-two genes lie within the support interval': the interval bounds are not printed")
def test_22_genes_in_support_interval():
    pass


def test_provenance_resolves(D, ds):
    codes = {c for c, d in zip(D["sources"]["code"], D["sources"]["ds"]) if d == ds}
    n = len(facts()) + 1
    for t in ("crosses", "qtl", "causal"):
        for d, pv in zip(D[t]["ds"], D[t]["pv"]):
            if d != ds:
                continue
            code, row = pv.split("#")
            assert code in codes and 2 <= int(row) <= n


def test_source_rows_kept(D, ds):
    ids = {r["id"] for r in rows(D["qtl"]) if r["ds"] == ds} | {r["id"] for r in rows(D["causal"]) if r["ds"] == ds}
    seen = set()
    for ch in R64:
        with open(os.path.join(ROOT, "docs", "v", f"d08_{ch}.json"), encoding="utf-8") as f:
            v = json.load(f)
        for rid, row in v["rows"]["d08:TF"].items():
            assert len(row) == len(v["cols"]["d08:TF"])
            seen.add(rid)
    assert seen == ids
    with open(os.path.join(ROOT, "docs", "v", "d08_extra.json"), encoding="utf-8") as f:
        x = json.load(f)
    assert len(x["d08:TF"]["rows"]) == 31 and len(x["d08:ART"]["rows"]) > 50
