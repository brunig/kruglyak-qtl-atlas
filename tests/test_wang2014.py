"""Validate the ingested Wang & Kruglyak 2014 (haloperidol) tables against the paper and its files.

    python3 -I -m pytest -q tests/test_wang2014.py      (after ingest/wang2014.py)

Claims from the article text are asserted as stated; where the supporting tables disagree, the
test is xfail(strict=True) with file and column in the reason.
"""
import collections
import importlib.util
import json
import os
import re
import html

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "sources", "wang2014")
EXT = os.path.join(SRC, "extracted")
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
    return D["datasets"]["id"].index("wang2014")


@pytest.fixture(scope="module")
def Q(D, ds):
    return [r for r in rows(D["qtl"]) if r["ds"] == ds]


@pytest.fixture(scope="module")
def T(D, ds):
    return {i: r for i, r in enumerate(rows(D["traits"])) if r["ds"] == ds}


def tsv(name):
    with open(os.path.join(EXT, name), encoding="utf-8") as f:
        return [l.rstrip("\n").split("\t") for l in f]


def s1_by_dose():
    out = collections.defaultdict(dict)
    for r in tsv("S1_Table.tsv")[1:]:
        if r[4] and r[5] != "df":
            out[r[3].split()[0]][r[4]] = r
    return out


def article_text():
    t = open(os.path.join(SRC, "article.html"), encoding="utf-8", errors="replace").read()
    t = re.sub(r"<(script|style)[^>]*>.*?</\1>", " ", t, flags=re.S)
    return re.sub(r"\s+", " ", html.unescape(re.sub(r"<[^>]+>", " ", t)))


# ---------------------------------------------------------------- sources ----------
def test_flattened_tables_match_the_docx_files():
    spec = importlib.util.spec_from_file_location("ext", os.path.join(EXT, "extract_tables.py"))
    ext = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ext)
    for name, f in ext.FILES.items():
        fresh = ext.flatten(os.path.join(SRC, "raw", f))
        on_disk = [r[: len(r)] for r in tsv(name + ".tsv")[1:]]
        assert [r + [""] * (len(on_disk[0]) - len(r)) for r in fresh] == on_disk, name


def test_text_facts_are_verbatim():
    t = article_text()
    for r in tsv("text_facts.tsv")[1:]:
        assert r[1] in t, r[1][:60]


# ---------------------------------------------------------------- design -----------
def test_one_cross_1008_segregants(D, ds):
    X = [r for r in rows(D["crosses"]) if r["ds"] == ds]
    assert len(X) == 1 and X[0]["n"] == 1008
    assert [D["parents"]["name"][X[0]["p1"]], D["parents"]["name"][X[0]["p2"]]] == ["BY", "RM"]


def test_five_doses(T):
    assert sorted(t["cond"] for t in T.values()) == sorted(["40 µM", "80 µM", "120 µM", "160 µM", "200 µM"])
    assert all(t["cat"] == "drug/chemical" for t in T.values())


def test_nine_distinct_loci(Q):
    labels = {r["k"].split()[1] for r in Q}
    assert len(labels) == 9


def test_qtl_rows_equal_S1_rows(Q):
    s1 = s1_by_dose()
    assert len(Q) == sum(len(v) - (1 if "Total" in v else 0) for v in s1.values()) == 32


def test_lod_above_calling_cutoff_in_the_additive_model(Q):
    # QTL were called at LOD 3 in the scans; S1 Table's LODs are drop-one values in the joint additive model, so a few
    # fall below 3 (15-2 at 80 uM, 13 at 160 uM). Recorded, not a contradiction of the calling rule.
    low = sorted(round(r["lod"], 3) for r in Q if r["lod"] < 3)
    assert low == [1.629, 2.313]


# ---------------------------------------------------------------- coordinates ------
def test_chromosome_only(Q):
    assert all(r["a"] is None and r["e"] is None and r["p"] is None and r["c"] in R64 for r in Q)


def test_swh1_snps_inside_yar042w_on_r64(D):
    g = D["genes"]["id"].index("YAR042W")
    a, e = D["genes"]["a"][g], D["genes"]["e"][g]
    pos = [int(r[4]) for r in tsv("S2_Table.tsv")[1:] if r[4].isdigit()]
    assert len(pos) == 22 and all(a <= p <= e for p in pos)


# ---------------------------------------------------------------- paper's numbers --
def test_additive_model_explains_35_9_to_54_7_percent():
    tot = {r[3].split()[0]: float(r[8]) for r in tsv("S1_Table.tsv")[1:] if r[4] == "Total"}
    four = [tot[d] for d in ("40", "80", "120", "160")]
    assert round(min(four), 1) == 35.9 and round(max(four), 1) == 54.7


@pytest.mark.xfail(strict=True, reason=(
    "Text: the chrI right-arm locus alone explained 21.4%, 26.7%, 21.4%, 11.6% at 40-160 uM. S1 Table "
    "(journal.pgen.1004894.s002.docx, row 1R, column %var) gives 21.3477, 26.6095, 21.390, 11.6005: the first two "
    "round to 21.3 and 26.6."))
def test_chrI_right_arm_variance_as_stated():
    s1 = s1_by_dose()
    assert [round(float(s1[d]["1R"][8]), 1) for d in ("40", "80", "120", "160")] == [21.4, 26.7, 21.4, 11.6]


def test_chrI_right_arm_variance_as_tabulated():
    s1 = s1_by_dose()
    assert [round(float(s1[d]["1R"][8]), 1) for d in ("40", "80", "120", "160")] == [21.3, 26.6, 21.4, 11.6]


def test_chrVII_and_chrXII_only_at_the_two_lower_doses():
    s1 = s1_by_dose()
    for lab in ("7", "12"):
        assert sorted(d for d in s1 if lab in s1[d]) == ["40", "80"]


@pytest.mark.xfail(strict=True, reason=(
    "Text: the chrXIII locus was only detected at the two lower doses. S1 Table (journal.pgen.1004894.s002.docx, "
    "160 uM table, row 13) also lists chrXIII at 160 uM (drop-one LOD 2.313, 0.69% variance)."))
def test_chrXIII_only_at_the_two_lower_doses():
    s1 = s1_by_dose()
    assert sorted(d for d in s1 if "13" in s1[d]) == ["40", "80"]


@pytest.mark.xfail(strict=True, reason=(
    "Text: 'RM alleles at loci on chromosomes V, XII and XV also confer greater resistance'. S1 Table "
    "(journal.pgen.1004894.s002.docx, first column, all doses) has no chromosome V locus; its nine loci are 1L, 1R, 2, "
    "7, 12, 13, 14, 15-1, 15-2."))
def test_a_chrV_locus_exists():
    s1 = s1_by_dose()
    assert any(lab.split("-")[0] == "5" for d in s1 for lab in s1[d])


def test_direction_reading_four_RM_five_BY(D, Q):
    up = {r["k"].split()[1]: D["parents"]["name"][r["up"]] for r in Q}
    assert sum(v == "RM" for v in up.values()) == 4 and sum(v == "BY" for v in up.values()) == 5
    assert up["1R"] == "RM" and up["14"] == "BY"


def test_interaction_counts_9_10_3_5_5():
    n = collections.Counter(r[3].split()[0] for r in tsv("S3_Table.tsv")[1:] if ":" in r[4])
    assert [n[d] for d in ("40", "80", "120", "160", "200")] == [9, 10, 3, 5, 5]


# ---------------------------------------------------------------- candidate genes --
def test_confirmed_genes_on_their_chromosomes(D, ds):
    C = [r for r in rows(D["causal"]) if r["ds"] == ds]
    assert {r["name"] for r in C} == {"SWH1", "MKT1", "IRA2"} and all(r["tier"] == "x" for r in C)
    q = D["qtl"]
    for r in C:
        assert q["c"][r["q1"]] == r["c"]


def test_swh1_on_the_right_arm(D):
    g = D["genes"]["n"].index("SWH1")
    assert D["genes"]["a"][g] > 151000   # chrI centromere CEN1 lies near 151.5 kb on R64


# ---------------------------------------------------------------- provenance -------
def test_provenance_resolves(D, ds):
    codes = {c for c, d in zip(D["sources"]["code"], D["sources"]["ds"]) if d == ds}
    for t in ("crosses", "qtl", "causal", "vexp"):
        for r in rows(D[t]):
            if r.get("ds", ds) != ds and t != "vexp":
                continue
            if t == "vexp" and D["traits"]["ds"][r["t"]] != ds:
                continue
            code, row = r["pv"].split("#")
            assert code in codes and int(row) >= 2


def test_original_rows_kept(D, ds):
    ids = {r["id"] for r in rows(D["qtl"]) if r["ds"] == ds}
    seen = set()
    for ch in R64:
        with open(os.path.join(ROOT, "docs", "v", f"w14_{ch}.json"), encoding="utf-8") as f:
            v = json.load(f)
        seen |= set(v["rows"]["w14:S1"])
        assert all(len(r) == len(v["cols"]["w14:S1"]) for r in v["rows"]["w14:S1"].values())
    assert seen == ids
