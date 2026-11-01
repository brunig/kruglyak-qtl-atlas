"""Validate the ingested Bloom et al. 2013 tables (Nature 494:234).

    python3 -I -m pytest -q tests/test_bloom2013.py      (after ./build.sh)

The article page is paywalled: only the abstract is public, so the one number checkable
against the paper's text is the 46 traits. Everything else checks the supplementary tables
against each other (Table S1's per-trait QTL counts against Table S3's rows), the genome
build (genes listed under each interval against R64-5-1 coordinates), coordinate sanity,
provenance, and that every source row is kept. Anomalies recorded in OPEN_QUESTIONS.md are
pinned so a rebuild that changes them is noticed.
"""
import collections
import json
import os

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "sources", "bloom2013", "raw")
S1 = os.path.join(RAW, "41586_2013_BFnature11867_MOESM86_ESM.xls")
S3 = os.path.join(RAW, "41586_2013_BFnature11867_MOESM88_ESM.xls")
R64 = {"chrI": 230218, "chrII": 813184, "chrIII": 316620, "chrIV": 1531933, "chrV": 576874, "chrVI": 270161,
       "chrVII": 1090940, "chrVIII": 562643, "chrIX": 439888, "chrX": 745751, "chrXI": 666816, "chrXII": 1078177,
       "chrXIII": 924431, "chrXIV": 784333, "chrXV": 1091291, "chrXVI": 948066}


def rows(tab):
    keys = list(tab)
    return [dict(zip(keys, v)) for v in zip(*(tab[k] for k in keys))]


def sheet(path):
    import xlrd
    sh = xlrd.open_workbook(path).sheet_by_index(0)
    return [sh.row_values(i) for i in range(1, sh.nrows)]


@pytest.fixture(scope="module")
def D():
    with open(os.path.join(ROOT, "docs", "data.json"), encoding="utf-8") as f:
        return json.load(f)


@pytest.fixture(scope="module")
def ds(D):
    return D["datasets"]["id"].index("bloom2013")


@pytest.fixture(scope="module")
def Q(D, ds):
    return [r for r in rows(D["qtl"]) if r["ds"] == ds]


@pytest.fixture(scope="module")
def T(D, ds):
    return [r for r in rows(D["traits"]) if r["ds"] == ds]


@pytest.fixture(scope="module")
def genes(D):
    return {g: (c, a, e) for c, a, e, g in zip(D["genes"]["c"], D["genes"]["a"], D["genes"]["e"], D["genes"]["id"])}


# ---------------------------------------------------------------- the paper -------
def test_46_traits(T):
    """Abstract: '46 quantitative traits'."""
    assert len(T) == 46


def test_one_cross_by_rm(D, ds):
    X = [r for r in rows(D["crosses"]) if r["ds"] == ds]
    assert len(X) == 1
    assert (D["parents"]["name"][X[0]["p1"]], D["parents"]["name"][X[0]["p2"]]) == ("BY", "RM")


def test_segregant_count_is_table_s1_maximum(D, ds):
    X = [r for r in rows(D["crosses"]) if r["ds"] == ds][0]
    assert X["n"] == max(int(r[3]) for r in sheet(S1)) == 1007


# ---------------------------------------------------------------- the tables ------
def test_591_qtl(Q):
    assert len(Q) == 591
    assert all(r["m"] == "w" for r in Q)


def test_qtl_per_trait_match_table_s1(D, Q):
    """Table S1 'Additve QTL detected' equals the number of Table S3 rows for every trait."""
    per = collections.Counter(D["traits"]["raw"][r["t"]] for r in Q)
    s1 = {r[0]: int(r[11]) for r in sheet(S1)}
    assert per == collections.Counter(s1)
    assert sum(s1.values()) == 591


def test_no_signed_effects(Q):
    assert all(r["ef"] is None and r["up"] == -1 for r in Q)
    assert all(0 < r["ve"] <= 1 for r in Q)


def test_variance_explained_table(D, ds):
    n = sum(1 for t in D["vexp"]["t"] if D["traits"]["ds"][t] == ds)
    assert n == 46


# ---------------------------------------------------------------- coordinates -----
def test_intervals_inside_chromosomes(Q):
    for r in Q:
        assert r["c"] in R64 and 1 <= r["a"] and r["e"] <= R64[r["c"]], r["id"]


def test_start_le_peak_le_end(Q):
    assert not [r["id"] for r in Q if not (r["a"] <= r["p"] <= r["e"])]


def test_build_r64_listed_genes_sit_under_their_interval(genes):
    """Every SGD systematic name Table S3 lists under an interval lies within 5 kb of that
    interval on R64-5-1 coordinates: the table's positions are R64-compatible."""
    roman = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI"]
    found = outside = 0
    for r in sheet(S3):
        ch = "chr" + roman[int(r[2]) - 1]
        for g in str(r[8]).split("|"):
            if g not in genes:
                continue
            c, a, e = genes[g]
            found += 1
            outside += not (c == ch and e >= r[5] - 5000 and a <= r[6] + 5000)
    assert found == 29205 and outside == 0


# ---------------------------------------------------------------- recorded anomalies
def test_some_lods_below_any_detection_threshold(Q):
    """OPEN_QUESTIONS (Bloom 2013, Q2): the LOD column is not a detection LOD."""
    lods = [r["lod"] for r in Q]
    assert sum(l < 2.5 for l in lods) == 34 and sum(l < 1 for l in lods) == 5 and min(lods) < 0.001


def test_sds_dose_as_displayed(T):
    """OPEN_QUESTIONS (Bloom 2013, Q4): Table S1 stores 0.1 with a '0%' format, displayed '10%'."""
    sds = [t for t in T if t["raw"] == "SDS"][0]
    assert sds["cond"] == "10%"


def test_condition_keys_link_to_bloom2019(D):
    if "bloom2019" not in D["datasets"]["id"]:
        pytest.skip("bloom2019 not built")
    by = collections.defaultdict(set)
    for ck, d in zip(D["traits"]["ck"], D["traits"]["ds"]):
        if ck:
            by[ck].add(D["datasets"]["id"][d])
    shared = {k for k, v in by.items() if {"bloom2013", "bloom2019"} <= v}
    assert {"6-azauracil 50ug/mL", "zeocin 25ug/mL", "neomycin 5mg/mL", "lithium chloride 100mM",
            "diamide 1.5mM", "formamide 2%", "ypd 37", "ynb ph8"} <= shared


# ---------------------------------------------------------------- provenance ------
def test_provenance_resolves(D, ds):
    codes = {c for c, d in zip(D["sources"]["code"], D["sources"]["ds"]) if d == ds}
    assert codes == {"b13:S1", "b13:S2", "b13:S3"}
    for t in ("crosses", "qtl", "vexp"):
        for pv, owner in zip(D[t]["pv"], range(len(D[t]["pv"]))):
            if pv and pv.startswith("b13:"):
                code, row = pv.split("#")
                assert code in codes and int(row) >= 2


def test_original_rows_kept_verbatim(Q):
    seen = {}
    for ch in R64:
        with open(os.path.join(ROOT, "docs", "v", f"b13_{ch}.json"), encoding="utf-8") as f:
            v = json.load(f)
        assert len(v["cols"]["b13:S3"]) == 9
        seen.update(v["rows"]["b13:S3"])
    assert set(seen) == {r["id"] for r in Q} and all(len(x) == 9 for x in seen.values())
    with open(os.path.join(ROOT, "docs", "v", "b13_extra.json"), encoding="utf-8") as f:
        ex = json.load(f)
    assert len(ex["b13:S1"]["rows"]) == 46 and len(ex["b13:S2"]["rows"]) == 17
