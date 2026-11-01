"""Validate the ingested Treusch et al. 2015 tables against the paper's text and against themselves.

    python3 -I -m pytest -q tests/test_treusch2015.py      (after python3 -I ingest/treusch2015.py)

Numbers quoted from the article (sources/treusch2015/article.html) are asserted as stated; where
the tables disagree the test is a strict xfail naming file and column. The data are not changed.
"""
import collections
import json
import os
import statistics

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "sources", "treusch2015", "raw")
R64 = {"chrI": 230218, "chrII": 813184, "chrIII": 316620, "chrIV": 1531933, "chrV": 576874, "chrVI": 270161,
       "chrVII": 1090940, "chrVIII": 562643, "chrIX": 439888, "chrX": 745751, "chrXI": 666816, "chrXII": 1078177,
       "chrXIII": 924431, "chrXIV": 784333, "chrXV": 1091291, "chrXVI": 948066}
S4 = "journal.pgen.1004913.s008.xlsx (S4 Table)"


def rows(tab):
    keys = list(tab)
    return [dict(zip(keys, v)) for v in zip(*(tab[k] for k in keys))]


@pytest.fixture(scope="module")
def D():
    with open(os.path.join(ROOT, "docs", "data.json"), encoding="utf-8") as f:
        return json.load(f)


@pytest.fixture(scope="module")
def ds(D):
    return D["datasets"]["id"].index("treusch2015")


@pytest.fixture(scope="module")
def Q(D, ds):
    return [r for r in rows(D["qtl"]) if r["ds"] == ds]


@pytest.fixture(scope="module")
def X(D, ds):
    return {i: r for i, r in enumerate(rows(D["crosses"])) if r["ds"] == ds}


@pytest.fixture(scope="module")
def grouped(Q):
    return [r for r in Q if r["k"] == "grouped round-robin QTL"]


def by_cond(D, L, raw):
    return [r for r in L if D["traits"]["raw"][r["t"]] == raw]


def rr_crosses(X):
    return [i for i, x in X.items() if x["code"] != "BYxRM"]


# ---------------------------------------------------------------- design ------------
def test_twelve_round_robin_crosses_each_strain_twice(X, D):
    rr = [X[i] for i in rr_crosses(X)]
    assert len(rr) == 12
    c = collections.Counter([x["p1"] for x in rr] + [x["p2"] for x in rr])
    assert len(c) == 12 and set(c.values()) == {2}


def test_65_strains_phenotyped_12_in_bold():
    import openpyxl
    ws = openpyxl.load_workbook(os.path.join(RAW, "journal.pgen.1004913.s005.xlsx")).active
    data = [r for r in ws.iter_rows(min_row=4) if r[0].value not in (None, "")]
    assert len(data) == 65
    assert sum(1 for r in data if r[0].font.bold) == 12


def test_build_is_sacCer3(D, ds):
    assert "sacCer3" in D["datasets"]["build"][ds]


# ---------------------------------------------------------------- grouped RR QTL ----
def test_60_regions_37_salt_23_caffeine(D, grouped):
    assert len(grouped) == 60
    assert len(by_cond(D, grouped, "NaCl")) == 37 and len(by_cond(D, grouped, "Caffeine")) == 23


def test_33_of_60_regions_in_one_cross(grouped):
    assert sum(1 for r in grouped if len(r["xs"]) == 1) == 33


def test_xs_matches_number_of_crosses_column(grouped):
    assert all(len(r["xs"]) == r["cc"] for r in grouped)


def _per_cross(D, X, L):
    return [sum(1 for r in L if x in r["xs"]) for x in rr_crosses(X)]


def test_salt_qtl_per_cross_mean_6_2_range_1_13(D, X, grouped):
    n = _per_cross(D, X, by_cond(D, grouped, "NaCl"))
    assert round(statistics.mean(n), 1) == 6.2 and (min(n), max(n)) == (1, 13)


def test_caffeine_qtl_per_cross_mean_3_7_range_0_7(D, X, grouped):
    n = _per_cross(D, X, by_cond(D, grouped, "Caffeine"))
    assert round(statistics.mean(n), 1) == 3.7 and (min(n), max(n)) == (0, 7)


@pytest.mark.xfail(strict=True, reason=(
    "Paper: 155 QTL identified across the twelve crosses (Results; Methods: '311 of 418 individually identified QTL ... "
    "resulting in 155 jointly identified QTL'). " + S4 + ", column '# Crosses with QTL', sums to 118 over the 60 regions. "
    "The 155 count per-cross QTL before doses of the same condition were combined; those per-cross rows are not in the bundle."))
def test_155_per_cross_qtl(grouped):
    assert sum(len(r["xs"]) for r in grouped) == 155


def test_per_cross_qtl_sum_as_observed(grouped):
    assert sum(len(r["xs"]) for r in grouped) == 118


@pytest.mark.xfail(strict=True, reason=(
    "Paper: 13 QTL regions recurred between salt and caffeine. Overlapping 2-LOD intervals in " + S4 +
    " (Chromosome, Lower/Upper bound of QTL) give 14 salt regions that overlap a caffeine region (and 14 the other way). "
    "The paper does not say how recurrence was defined."))
def test_13_regions_shared_between_conditions(D, grouped):
    sal, caf = by_cond(D, grouped, "NaCl"), by_cond(D, grouped, "Caffeine")
    assert sum(any(a["c"] == b["c"] and a["a"] <= b["e"] and b["a"] <= a["e"] for b in caf) for a in sal) == 13


def _with_gene(D, grouped, raw, gene):
    g = D["genes"]["n"].index(gene)
    return [r for r in by_cond(D, grouped, raw) if r["g"] == g]


def test_ena_in_11_of_12_crosses(D, grouped):
    r = [q for q in by_cond(D, grouped, "NaCl") if q["c"] == "chrIV" and q["a"] <= 509600 <= q["e"]]
    assert len(r) == 1 and len(r[0]["xs"]) == 11


def test_tor1_caffeine_in_8_crosses(D, grouped):
    r = _with_gene(D, grouped, "Caffeine", "TOR1")
    assert len(r) == 1 and len(r[0]["xs"]) == 8 and r[0]["gt"] == "x"


def test_tor2_caffeine_5_crosses_salt_3_crosses(D, grouped):
    assert [len(r["xs"]) for r in _with_gene(D, grouped, "Caffeine", "TOR2")] == [5]
    assert [len(r["xs"]) for r in _with_gene(D, grouped, "NaCl", "TOR2")] == [3]


def test_whi2_confirmed_in_crosses_7_and_8(D, X, grouped):
    code = {i: x["code"] for i, x in X.items()}
    for raw in ("NaCl", "Caffeine"):
        r = _with_gene(D, grouped, raw, "WHI2")
        assert len(r) == 1 and r[0]["gt"] == "x" and {"7", "8"} <= {code[x] for x in r[0]["xs"]}


def test_mss4_caffeine_in_crosses_2_and_3(D, X, grouped):
    code = {i: x["code"] for i, x in X.items()}
    r = _with_gene(D, grouped, "Caffeine", "MSS4")
    assert len(r) == 1 and sorted(code[x] for x in r[0]["xs"]) == ["2", "3"]


def test_every_grouped_qtl_has_allele_strains(grouped):
    assert all(r["ups"] and r["dns"] for r in grouped)


def test_allele_strains_are_parents_of_detecting_crosses(D, X, grouped):
    """S7 strains should be parents of the crosses in which the QTL was detected."""
    bad = []
    for r in grouped:
        parents = {p for x in r["xs"] for p in (X[x]["p1"], X[x]["p2"])}
        if not set(r["ups"]) | set(r["dns"]) <= parents:
            bad.append(r["id"])
    assert not bad, bad


# ---------------------------------------------------------------- MAT QTL --------------
def _mat(D, Q):
    return [r for r in Q if r["k"].startswith("MAT-linked")]


def test_mat_qtl_6_in_ypd_19_under_stress(D, Q):
    m = _mat(D, Q)
    ypd = [r for r in m if D["traits"]["raw"][r["t"]] == "MAT:YPD"]
    assert len(ypd) == 6 and len(m) - len(ypd) == 19


def test_14_of_19_stress_mat_qtl_in_multiple_conditions(D, Q):
    st = [r for r in _mat(D, Q) if D["traits"]["raw"][r["t"]] != "MAT:YPD"]
    multi = sum(any(o is not r and o["x"] == r["x"] and o["c"] == r["c"] and o["a"] <= r["e"] and r["a"] <= o["e"] for o in st) for r in st)
    assert multi == 14


def test_gpa1_mat_qtl_in_two_crosses(D, X, Q):
    g = D["genes"]["n"].index("GPA1")
    ypd = [r for r in _mat(D, Q) if r["g"] == g and D["traits"]["raw"][r["t"]] == "MAT:YPD"]
    assert sorted(X[r["x"]]["code"] for r in ypd) == ["10", "5"]


# ---------------------------------------------------------------- coordinates ---------
def test_intervals_inside_chromosomes(Q):
    for r in Q:
        assert r["c"] in R64 and 1 <= r["a"] and r["e"] <= R64[r["c"]], r["id"]


def test_start_le_peak_le_end(Q):
    assert not [r["id"] for r in Q if not (r["a"] <= r["p"] <= r["e"])]


def test_lod_threshold_5(Q):
    assert min(r["lod"] for r in Q) > 5


def _outside(D, Q):
    G = D["genes"]
    return [r["id"] for r in Q if r["g"] >= 0 and not (G["c"][r["g"]] == r["c"] and G["a"][r["g"]] <= r["e"] and G["e"][r["g"]] >= r["a"])]


@pytest.mark.xfail(strict=True, reason=(
    "Two candidate genes named by the paper lie outside the 2-LOD interval of their row (SGD R64-5-1 bounds): "
    "CYR1 (chrX:425,157-431,237) for S4 Caff.chr10.2 (431,300-477,200), 63 bp short; and GPA1 "
    "(chrVIII:113,375-115,043) for S8 Cross 10, 15mM Caffeine, chr08 (83,500-103,900), 9.5 kb past the bound. "
    "Sources: journal.pgen.1004913.s008.xlsx and s012.xlsx, columns Lower/Upper bound of QTL and Candidate QTG."))
def test_candidate_gene_inside_interval(D, Q):
    assert not _outside(D, Q)


def test_candidate_genes_outside_interval_are_exactly_two(D, Q):
    assert sorted(_outside(D, Q)) == ["t15-x-S4-12", "t15-x-S8-24"]


# ---------------------------------------------------------------- provenance ----------
def test_provenance_resolves(D, ds):
    codes = {c for c, d in zip(D["sources"]["code"], D["sources"]["ds"]) if d == ds}
    for t in ("crosses", "qtl", "causal"):
        for pv, d in zip(D[t]["pv"], D[t]["ds"]):
            if d == ds:
                code, row = pv.split("#")
                assert code in codes and int(row) >= 2


def test_source_codes_globally_unique(D):
    assert len(D["sources"]["code"]) == len(set(D["sources"]["code"]))


def test_original_rows_kept(D, ds, Q):
    seen = collections.Counter()
    for ch in R64:
        with open(os.path.join(ROOT, "docs", "v", f"t15_{ch}.json"), encoding="utf-8") as f:
            v = json.load(f)
        for code, rs in v["rows"].items():
            seen[code] += len(rs)
    assert seen == {"t15:S4": 60, "t15:S3": 12, "t15:S8": 25, "t15:S7": 60}
    with open(os.path.join(ROOT, "docs", "v", "t15_extra.json"), encoding="utf-8") as f:
        ex = json.load(f)
    assert sorted(ex) == sorted(f"t15:S{i}" for i in range(1, 11))
