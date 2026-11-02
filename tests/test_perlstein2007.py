"""Validate the ingested Perlstein et al. 2007 tables against the paper (abstract and supplement).

    python3 -I -m pytest -q tests/test_perlstein2007.py      (after the ingest)

Paper claims are asserted as stated; where the files disagree the test is a strict xfail naming
file and column. The main text is paywalled, so the checkable claims are the abstract's and the
supplementary legends'.
"""
import csv
import html
import json
import os
import re
import statistics

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "sources", "perlstein2007")
EXT = os.path.join(SRC, "extracted")
R64 = {"chrI": 230218, "chrII": 813184, "chrIII": 316620, "chrIV": 1531933, "chrV": 576874, "chrVI": 270161,
       "chrVII": 1090940, "chrVIII": 562643, "chrIX": 439888, "chrX": 745751, "chrXI": 666816, "chrXII": 1078177,
       "chrXIII": 924431, "chrXIV": 784333, "chrXV": 1091291, "chrXVI": 948066}


def rows(tab):
    keys = list(tab)
    return [dict(zip(keys, v)) for v in zip(*(tab[k] for k in keys))]


def ext(name):
    with open(os.path.join(EXT, name), encoding="utf-8") as f:
        return list(csv.DictReader((l for l in f if not l.startswith("#")), delimiter="\t"))


@pytest.fixture(scope="module")
def D():
    with open(os.path.join(ROOT, "docs", "data.json"), encoding="utf-8") as f:
        return json.load(f)


@pytest.fixture(scope="module")
def ds(D):
    return D["datasets"]["id"].index("perlstein2007")


@pytest.fixture(scope="module")
def Q(D, ds):
    return [r for r in rows(D["qtl"]) if r["ds"] == ds]


@pytest.fixture(scope="module")
def L(Q):
    return [r for r in Q if r["id"].startswith("p07-w-")]


@pytest.fixture(scope="module")
def H(Q):
    return [r for r in Q if r["id"].startswith("p07-h-")]


@pytest.fixture(scope="module")
def genes(D):
    return rows(D["genes"])


def norm(s):
    return re.sub(r"\s+", " ", s).strip()


def test_quotes_verbatim():
    import pdfplumber
    with open(os.path.join(SRC, "article.html"), encoding="utf-8", errors="replace") as f:
        art = norm(html.unescape(re.sub(r"<[^>]+>", " ", f.read())))
    pdf = {}
    for r in ext("text_facts.tsv"):
        if r["source"] == "article.html":
            src = art
        else:
            if r["source"] not in pdf:
                with pdfplumber.open(os.path.join(SRC, r["source"])) as p:
                    pdf[r["source"]] = norm(" ".join((pg.extract_text() or "") for pg in p.pages))
            src = pdf[r["source"]]
        assert r["quote"] in src, r["quote"]


def test_flattened_tables_match_legend_row_counts():
    # Table 3: rows 2-220; Table 4: rows 2-124; Table 5: rows 4-127 (legends)
    assert (len(ext("table_s3.tsv")), len(ext("table_s4.tsv")), len(ext("table_s5.tsv"))) == (219, 123, 124)
    with open(os.path.join(EXT, "unparsed_lines.tsv")) as f:
        assert len(f.read().strip().split("\n")) == 1   # header only: nothing left unparsed


def test_cross_and_segregants(D, ds):
    X = [r for r in rows(D["crosses"]) if r["ds"] == ds]
    assert len(X) == 1 and X[0]["n"] == 104
    P = [r for r in rows(D["parents"]) if r["ds"] == ds]
    assert {p["name"] for p in P} == {"BY", "RM"}


def test_124_distinct_linkages(L):
    assert len(L) == 124


def test_83_compounds_with_linkages():
    names = {r["compound"] for r in ext("table_s4.tsv")}
    assert len(names) == 84
    assert len({n.replace(" (CSM media)", "") for n in names}) == 83   # rapamycin in two media counted once


@pytest.mark.xfail(strict=True, reason=(
    "Abstract: 'a collection of 100 diverse small molecules'. Supplementary Table 1 "
    "(41588_2007_BFng1991_MOESM236_ESM.pdf, column 'SMP name') lists 92 numbered SMPs."))
def test_100_smps_listed():
    assert len(ext("table_s1.tsv")) == 100


def test_eight_hotspots(H):
    assert len(H) == 8
    assert sorted((h["c"], h["p"]) for h in H) == sorted([("chrI", 180000), ("chrIII", 80000), ("chrXII", 660000), ("chrXIII", 20000),
                                                           ("chrXIV", 440000), ("chrXIV", 480000), ("chrXV", 160000), ("chrXV", 240000)])


def test_hotspot_smp_lists_match_counts_except_chrXIII():
    hs = {r["column"]: int(r["n_smp"]) for r in ext("table_s6_hotspots.tsv")}
    got = {}
    for r in ext("table_s6_hotspot_smps.tsv"):
        got[r["column"]] = got.get(r["column"], 0) + 1
    assert {k: v for k, v in got.items() if k != "4"} == {k: v for k, v in hs.items() if k != "4"}


@pytest.mark.xfail(strict=True, reason=(
    "Supplementary Table 6 (41588_2007_BFng1991_MOESM241_ESM.pdf) gives 26 linked SMPs for the chrXIII 20000 hotspot "
    "(column D header) but its column lists 25 names (extracted/table_s6_hotspot_smps.tsv, column 4)."))
def test_chrXIII_hotspot_lists_26_smps():
    assert sum(1 for r in ext("table_s6_hotspot_smps.tsv") if r["column"] == "4") == 26


def test_coordinates_inside_R64(Q):
    for r in Q:
        assert r["c"] in R64 and 0 <= r["a"] <= r["p"] <= r["e"] <= R64[r["c"]], r["id"]


def test_build_consistent_with_R64(L, genes):
    by_id = {g["id"]: g for g in genes}
    d = []
    for r in L:
        orf = r["lv"].split("_")[0]
        orf = orf[:-2] if orf.endswith(".1") else orf
        if orf in by_id:
            g = by_id[orf]
            assert g["c"] == r["c"], r["id"]
            d.append(abs(r["p"] - (g["a"] + g["e"]) / 2))
    assert len(d) >= 90 and statistics.median(d) < 1000 and max(d) < 5000


def test_table5_labels_agree_with_their_time_points(D, ds):
    extra = json.load(open(os.path.join(ROOT, "docs", "v", "p07_extra.json"), encoding="utf-8"))["p07:ALIGN"]
    rs = [dict(zip(extra["cols"], r)) for r in extra["rows"]]
    assert len(rs) == 124
    for r in rs:
        m = re.search(r"AT\.(\d+)$", r["table5_trait"])
        if m:
            assert m.group(1) == r["time_h"], r
    # three Table 5 rows have no Table 4 row at the same time point: cinnarizine 44 h, H2O2 132 h, tyrothricin 80 h
    assert sorted(r["table5_trait"] for r in rs if r["label_source"] == "p07:S3") == ["H2O2HighAT.132", "cinnarazineAT.44", "tyrothricinAT.80"]


def test_candidates(D, ds, L, H, genes):
    C = [r for r in rows(D["causal"]) if r["ds"] == ds]
    assert sorted((c["name"], c["tier"]) for c in C) == sorted([("LEU2", "x"), ("PHO84", "x"), ("PHO84", "x"), ("SWH1", "c"), ("MKS1", "c"), ("MKS1", "c")])
    q = rows(D["qtl"])
    for c in C:
        qq = q[c["q1"]]
        assert qq["c"] == c["c"]
        if c["name"] == "LEU2":   # inside the niguldipine interval
            assert qq["a"] <= c["e"] and qq["e"] >= c["a"]


def test_PHO84_sits_just_outside_the_phenol_intervals(D, ds):
    """Recorded, not hidden: PHO84 (chrXIII 23,886-25,919) ends ~1.7 kb left of both intervals, which start at the peak marker."""
    q = rows(D["qtl"])
    for c in rows(D["causal"]):
        if c["ds"] == ds and c["name"] == "PHO84":
            qq = q[c["q1"]]
            assert qq["a"] == qq["p"] == 27644 and 0 < qq["a"] - c["e"] < 2000


def test_no_effects_invented(Q):
    assert all(r["ef"] is None and r["up"] == -1 for r in Q)


def test_provenance(D, ds):
    codes = set(D["sources"]["code"])
    for t in ("crosses", "qtl", "causal"):
        for r in rows(D[t]):
            if r.get("ds", ds) == ds and str(r["pv"]).startswith("p07:"):
                code, row = r["pv"].split("#")
                assert code in codes and int(row) >= 2


def test_verbatim_rows_kept(Q):
    seen = set()
    for c in R64:
        v = json.load(open(os.path.join(ROOT, "docs", "v", f"p07_{c}.json"), encoding="utf-8"))
        for code, rs in v["rows"].items():
            for rid, val in rs.items():
                assert len(val) == len(v["cols"][code])
                seen.add(rid)
    assert {r["id"] for r in Q} <= seen
    extra = json.load(open(os.path.join(ROOT, "docs", "v", "p07_extra.json"), encoding="utf-8"))
    assert len(extra["p07:S3"]["rows"]) == 219 and len(extra["p07:S1"]["rows"]) == 92
