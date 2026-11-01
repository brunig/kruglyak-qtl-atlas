"""Validate the ingested Bloom et al. 2015 tables against the paper and against themselves.

    python3 -I -m pytest -q tests/test_bloom2015.py      (after ./build.sh or python3 -I ingest/bloom2015.py)

Claims from the article text are asserted as stated; where the files disagree the test is a strict
xfail naming the file and column. Nothing is changed to match.
"""
import collections
import json
import os
import statistics

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "sources", "bloom2015")
RAW = os.path.join(SRC, "raw")
SD2 = "41467_2015_BFncomms9712_MOESM726_ESM.xls"
SD4 = "41467_2015_BFncomms9712_MOESM728_ESM.xls"
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
    return D["datasets"]["id"].index("bloom2015")


@pytest.fixture(scope="module")
def Q(D, ds):
    return [r for r in rows(D["qtl"]) if r["ds"] == ds]


def sheet(name, sh):
    import xlrd
    s = xlrd.open_workbook(os.path.join(RAW, name)).sheet_by_name(sh) if isinstance(sh, str) else \
        xlrd.open_workbook(os.path.join(RAW, name)).sheet_by_index(sh)
    return [s.row_values(i) for i in range(s.nrows)]


def interactions(sh):
    return [r for r in sheet(SD4, sh)[1:] if r[0] != ""]


def vc_column(col):
    """Supplementary Data 2, model y=Za+Zi+Zp+e (columns 12-17), rows of the 20 traits."""
    rs = sheet(SD2, 0)
    return [r[col] for r in rs[2:22]]


# ------------------------------------------------------------------ design ------------------
def test_one_cross_4390_segregants(D, ds):
    X = [r for r in rows(D["crosses"]) if r["ds"] == ds]
    assert len(X) == 1 and X[0]["n"] == 4390
    P = D["parents"]
    assert {P["name"][X[0]["p1"]], P["name"][X[0]["p2"]]} == {"BY", "RM"}


def test_20_traits(D, ds):
    assert sum(1 for d in D["traits"]["ds"] if d == ds) == 20


def _markers():
    with open(os.path.join(SRC, "extracted", "markers.tsv")) as f:
        return [l for l in f if not l.startswith("#")][1:]


def test_28220_markers_methods():
    # Methods: "28,220 unique single nucleotide polymorphisms were retained"
    assert len(_markers()) == 28220


@pytest.mark.xfail(strict=True, reason="Introduction says 28,820 variant sites and the Supplementary Data 5 caption says a 4390x28880 matrix; "
                   "the genotype matrix in Supplementary Data 5 (genotypes.RData, gdata) has 28,220 columns, as the Methods say.")
def test_28820_markers_introduction():
    assert len(_markers()) == 28820


def test_markers_inside_R64():
    for l in _markers():
        _, _, chrom, pos, _, _ = l.rstrip("\n").split("\t")
        assert 1 <= int(pos) <= R64[chrom]


# ------------------------------------------------------------------ additive QTL -------------
def test_797_additive_qtl(Q):
    assert len(Q) == 797


def test_per_trait_median_42_5_range_17_56(Q):
    n = collections.Counter(r["t"] for r in Q)
    assert len(n) == 20 and statistics.median(n.values()) == 42.5
    assert min(n.values()) == 17 and max(n.values()) == 56


def test_180_qtl_explain_1pct_median_ci_10_3kb(Q):
    big = [r for r in Q if r["ve"] >= 0.01]
    assert len(big) == 180
    assert round(statistics.median(r["e"] - r["a"] for r in big) / 1000, 1) == 10.3


def test_46_qtl_above_5pct_largest_26pct(Q):
    assert sum(1 for r in Q if r["ve"] > 0.05) == 46
    assert round(100 * max(r["ve"] for r in Q)) == 26


def test_median_effect_size_0_38pct(Q):
    assert round(100 * statistics.median(r["ve"] for r in Q), 2) == 0.38


@pytest.mark.xfail(strict=True, reason="Text: 'In 31 cases, QTL could be refined to a single gene (Supplementary Data 3)'. "
                   "Supplementary Data 3, column 'systematic name', lists exactly one gene for 25 QTL (13 distinct genes).")
def test_31_single_gene_qtl(D, ds):
    assert sum(1 for r in rows(D["causal"]) if r["ds"] == ds) == 31


def test_single_gene_qtl_as_observed(D, ds):
    C = [r for r in rows(D["causal"]) if r["ds"] == ds]
    assert len(C) == 25 and len({r["g"] for r in C}) == 13 and {r["tier"] for r in C} == {"c"}


def test_single_gene_inside_its_interval(D, ds):
    G = D["genes"]
    for r in rows(D["causal"]):
        if r["ds"] != ds:
            continue
        q = {k: D["qtl"][k][r["q1"]] for k in ("c", "a", "e")}
        assert G["c"][r["g"]] == q["c"] and G["a"][r["g"]] <= q["e"] and G["e"][r["g"]] >= q["a"], r["id"]


def test_effect_sign_and_direction(Q, D):
    P = D["parents"]["name"]
    for r in Q:
        assert P[r["up"]] == ("RM" if r["ef"] > 0 else "BY")


# ------------------------------------------------------------------ coordinates --------------
def test_inside_R64_and_a_le_p_le_e(Q):
    for r in Q:
        assert r["c"] in R64 and 1 <= r["a"] <= r["p"] <= r["e"] <= R64[r["c"]], r["id"]


def test_lead_marker_matches_source_position(Q):
    import xlrd
    s = xlrd.open_workbook(os.path.join(RAW, "41467_2015_BFncomms9712_MOESM727_ESM.xls")).sheet_by_index(0)
    for r in Q:
        row = int(r["pv"].split("#")[1])
        assert r["lv"].split("_")[:2] == [s.cell_value(row - 1, 0), str(int(s.cell_value(row - 1, 3)))]


def test_forward_peak_outside_interval_count(Q):
    # documented in datasets.interval: the forward-search marker lies outside the refineqtl interval for 11 QTL
    assert sum(1 for r in Q if not (r["a"] <= int(r["lv"].split("_")[1]) <= r["e"])) == 11


# ------------------------------------------------------------------ interactions --------------
def test_205_full_266_marginal():
    assert len(interactions("Full Scan")) == 205 and len(interactions("Marginal Scan")) == 266


def test_full_scan_153_36_16():
    c = collections.Counter(int(r[12]) for r in interactions("Full Scan"))
    assert (c[0], c[1], c[2]) == (153, 36, 16)


def test_largest_interaction_3_3pct_median_0_31pct():
    f = [r[9] for r in interactions("Full Scan")]
    assert round(100 * max(f), 1) == 3.3 and round(100 * statistics.median(f), 2) == 0.31


def test_about_13_interactions_per_trait():
    assert round(len(interactions("Marginal Scan")) / 20) == 13


def test_interaction_marker_indices_match_positions():
    mk = {}
    for l in _markers():
        i, _, chrom, pos, _, _ = l.rstrip("\n").split("\t")
        mk[int(i)] = (chrom, int(pos))
    for sh in ("Full Scan", "Marginal Scan"):
        for r in interactions(sh):
            assert mk[int(r[1])] == (r[3], int(r[4])) and mk[int(r[2])] == (r[5], int(r[6]))


def _both(d):
    F, M = interactions("Full Scan"), interactions("Marginal Scan")

    def near(a, b):
        return a[0] == b[0] and abs(a[1] - b[1]) <= d
    n = 0
    for f in F:
        a = ((f[3], f[4]), (f[5], f[6]))
        n += any(f[0] == m[0] and ((near(a[0], (m[3], m[4])) and near(a[1], (m[5], m[6]))) or
                                   (near(a[0], (m[5], m[6])) and near(a[1], (m[3], m[4])))) for m in M)
    return n


@pytest.mark.xfail(strict=True, reason="Text: '172 interactions detected by both approaches'. Supplementary Data 4 gives no matching rule; "
                   "same trait with both loci within 50 kb matches 170 full-scan rows, within 100 kb 173, identical markers 13, "
                   "identical closest additive QTL pair 165.")
def test_172_detected_by_both():
    assert _both(50000) == 172


# ------------------------------------------------------------------ variance components (Supplementary Data 2)
def test_additive_variance_8_6_to_70_4_median_43_3():
    v = vc_column(11)
    assert (round(100 * min(v), 1), round(100 * max(v), 1), round(100 * statistics.median(v), 1)) == (8.6, 70.4, 43.3)


def test_additive_mean_43_7():
    assert round(100 * statistics.mean(vc_column(11)), 1) == 43.7


@pytest.mark.xfail(strict=True, reason="Text: interaction variance 'ranged from 2.2 to 21.2%'. Supplementary Data 2, column 'Interactions (GxG)' "
                   "of the y=Za+Zi+Zp+e model: 2.29% to 21.26% (2.3 to 21.3 rounded; the text's values match truncation). Median 9.2% agrees.")
def test_interaction_variance_2_2_to_21_2():
    v = vc_column(13)
    assert (round(100 * min(v), 1), round(100 * max(v), 1)) == (2.2, 21.2)


def test_interaction_median_9_2():
    assert round(100 * statistics.median(vc_column(13)), 1) == 9.2


@pytest.mark.xfail(strict=True, reason="Text: residual strain repeatability 'ranged from 0.05 to 21.4%'. Supplementary Data 2, column "
                   "'Repeatibility' of the y=Za+Zi+Zp+e model: minimum 0.56%; maximum 21.4% and median 8.8% agree.")
def test_repeatability_0_05_to_21_4():
    v = vc_column(15)
    assert round(100 * min(v), 2) == 0.05 and round(100 * max(v), 1) == 21.4


def test_repeatability_max_and_median():
    v = vc_column(15)
    assert round(100 * max(v), 1) == 21.4 and round(100 * statistics.median(v), 1) == 8.8


# ------------------------------------------------------------------ provenance ----------------
def test_provenance_resolves_and_codes_unique(D, ds):
    S = D["sources"]
    assert len(S["code"]) == len(set(S["code"]))
    codes = {c for c, d in zip(S["code"], S["ds"]) if d == ds}
    for t in ("crosses", "qtl", "causal"):
        for r in rows(D[t]):
            if r.get("ds", ds) == ds and r["pv"].startswith("b15:"):
                assert r["pv"].split("#")[0] in codes


def test_original_rows_kept_verbatim(D):
    seen = collections.Counter()
    for ch in R64:
        with open(os.path.join(ROOT, "docs", "v", f"b15_{ch}.json"), encoding="utf-8") as f:
            v = json.load(f)
        w = len(v["cols"]["b15:SD3"])
        for rid, val in v["rows"]["b15:SD3"].items():
            assert len(val) == w
            seen["c" if "-c-" in rid else "w"] += 1
    assert seen == {"w": 797, "c": 25}
    with open(os.path.join(ROOT, "docs", "v", "b15_extra.json"), encoding="utf-8") as f:
        x = json.load(f)
    assert len(x["b15:SD4f"]["rows"]) == 206 and len(x["b15:SD4m"]["rows"]) == 267
