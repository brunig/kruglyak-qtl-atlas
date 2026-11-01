"""Validate the ingested Ehrenreich et al. 2012 tables against the paper and against themselves.

    python3 -I -m pytest -q tests/test_ehrenreich2012.py      (after python3 -I ingest/ehrenreich2012.py)

Numbers come from the article text (sources/ehrenreich2012/article.html) unless marked as Table S1,
which is a picture inside journal.pgen.1002570.s015.doc; its legible totals were transcribed by eye
and are labelled so. Disagreements are strict xfails naming the file and column.
"""
import collections
import csv
import itertools
import json
import os

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "sources", "ehrenreich2012", "raw")
S10 = os.path.join(RAW, "journal.pgen.1002570.s013.txt")
S11 = os.path.join(RAW, "journal.pgen.1002570.s014.txt")
P4 = ["BY", "RM", "YJM789", "YPS163"]


def rows(tab):
    keys = list(tab)
    return [dict(zip(keys, v)) for v in zip(*(tab[k] for k in keys))]


@pytest.fixture(scope="module")
def D():
    with open(os.path.join(ROOT, "docs", "data.json"), encoding="utf-8") as f:
        return json.load(f)


@pytest.fixture(scope="module")
def ds(D):
    return D["datasets"]["id"].index("ehrenreich2012")


@pytest.fixture(scope="module")
def Q(D, ds):
    return [r for r in rows(D["qtl"]) if r["ds"] == ds]


@pytest.fixture(scope="module")
def X(D, ds):
    return [r for r in rows(D["crosses"]) if r["ds"] == ds]


@pytest.fixture(scope="module")
def CG(D, ds):
    return [r for r in rows(D["causal"]) if r["ds"] == ds]


def pname(D, i):
    return D["parents"]["name"][i]


def tname(D, r):
    return D["traits"]["raw"][r["t"]]


# ---------------------------------------------------------------- design ----------
def test_four_parents(D, ds):
    assert sorted(n for n, d in zip(D["parents"]["name"], D["parents"]["ds"]) if d == ds) == sorted(P4)


def test_all_six_pairwise_crosses(D, X):
    pairs = {frozenset((pname(D, x["p1"]), pname(D, x["p2"]))) for x in X}
    assert len(X) == 6 and pairs == {frozenset(p) for p in itertools.combinations(P4, 2)}


def test_13_chemicals_78_experiments(D, ds):
    T = [r for r in rows(D["traits"]) if r["ds"] == ds]
    assert len(T) == 13
    # every trait has a dose for each of the 6 crosses (File S1): 78 X-QTL experiments
    assert all(t["cond"].count(";") == 5 for t in T) and 13 * 6 == 78


# ---------------------------------------------------------------- loci ------------
def test_411_distinct_loci(Q):
    assert len(Q) == 411 and all(q["m"] == "x" for q in Q)


def test_32_loci_per_condition_min_8_cycloheximide_max_57_zeocin(D, Q):
    n = collections.Counter(tname(D, q) for q in Q)
    assert round(len(Q) / len(n)) == 32
    assert min(n.values()) == n["Cycloheximide"] == 8
    assert max(n.values()) == n["Zeocin"] == 57


@pytest.mark.xfail(strict=True, reason=(
    "Text: 837 peaks across the 78 X-QTL experiments. File S10 (journal.pgen.1002570.s013.txt), columns "
    "BYxRM ... YJMxYPS, has 836 TRUE values (one per detected per-cross peak)."))
def test_837_peaks(Q):
    assert sum(len(q["xs"]) for q in Q) == 837


def test_peaks_as_observed_and_mean_per_experiment(Q):
    n = sum(len(q["xs"]) for q in Q)
    assert n == 836
    assert round(n / 78, 1) == 10.7   # the text's "average of 10.7 peaks per trait per cross" holds with 836 too


def test_copper_byxrm_18_by_2_rm(D, Q):
    cu = [q for q in Q if tname(D, q) == "Copper sulfate"]
    byxrm = D["crosses"]["id"].index("e12-BYxRM")
    fav = collections.Counter(pname(D, q["ups"][q["xs"].index(byxrm)]) for q in cu if byxrm in q["xs"])
    assert fav == {"BY": 18, "RM": 2}


def test_ups_dns_pair_with_xs(D, Q):
    for q in Q:
        assert len(q["xs"]) == len(q["ups"]) == len(q["dns"]) >= 1
        for x, u, d in zip(q["xs"], q["ups"], q["dns"]):
            assert {u, d} == {D["crosses"]["p1"][x], D["crosses"]["p2"][x]}


# ------------------------------------------------- allelic singletons and doubletons
def _classes():
    """Classify File S10 loci with the Table S1 rules: exact patterns or one missing peak."""
    X6 = ["BYxRM", "BYxYJM", "BYxYPS", "RMxYJM", "RMxYPS", "YJMxYPS"]
    P = ["BY", "RM", "YJM", "YPS"]
    out = []
    with open(S10, newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            fav = {}
            for c in X6:
                s = int(r[c + "_sign"]); a, b = c.split("x")
                fav[c] = None if s == 0 else a if s == 1 else b
            cl = set()
            for S in P:
                inv = [c for c in X6 if S in c.split("x")]
                if any(fav[c] for c in X6 if c not in inv):
                    continue
                det = [c for c in inv if fav[c]]
                if len(det) >= 2 and all(fav[c] == S for c in det):
                    cl.add((S, "singleton", "resistant"))
                elif len(det) >= 2 and all(fav[c] != S for c in det):
                    cl.add((S, "singleton", "susceptible"))
            for A, B in [("BY", "RM"), ("BY", "YJM"), ("BY", "YPS")]:   # each doubleton split once
                C, E = [p for p in P if p not in (A, B)]
                inv = [c for c in X6 if set(c.split("x")) in ({A, C}, {A, E}, {B, C}, {B, E})]
                if any(fav[c] for c in X6 if c not in inv):
                    continue
                det = [c for c in inv if fav[c]]
                if len(det) >= 3 and all(fav[c] in (A, B) for c in det):
                    cl.add((A + "=" + B, "doubleton", "first pair resistant"))
                elif len(det) >= 3 and all(fav[c] in (C, E) for c in det):
                    cl.add((A + "=" + B, "doubleton", "second pair resistant"))
            out.append(cl)
    return out


def test_135_singletons_28_doubletons():
    cl = _classes()
    assert sum(any(k[1] == "singleton" for k in c) for c in cl) == 135
    assert sum(any(k[1] == "doubleton" for k in c) for c in cl) == 28
    assert all(len(c) <= 1 for c in cl)


def test_remaining_248_loci_are_210_unclassified_plus_38_series():
    # text: 29 series with at least 3 alleles, 9 needing 4 alleles, and 210 loci (51%) not interpretable
    cl = _classes()
    other = sum(not c for c in cl)
    assert other == 248 == 210 + 29 + 9
    assert round(100 * 210 / 411) == 51


def test_table_s1_class_totals_transcribed():
    """Table S1 (image in journal.pgen.1002570.s015.doc), legible totals transcribed by eye."""
    n = collections.Counter(k for c in _classes() for k in c)
    tot = lambda s: n[(s, "singleton", "resistant")] + n[(s, "singleton", "susceptible")]
    assert (n[("BY", "singleton", "resistant")], n[("BY", "singleton", "susceptible")], tot("BY")) == (26, 16, 42)
    assert (n[("YJM", "singleton", "resistant")], n[("YJM", "singleton", "susceptible")], tot("YJM")) == (7, 11, 18)
    assert (n[("YPS", "singleton", "resistant")], n[("YPS", "singleton", "susceptible")], tot("YPS")) == (9, 17, 26)
    assert tot("RM") == 135 - 42 - 18 - 26 == 49
    dt = lambda p: n[(p, "doubleton", "first pair resistant")] + n[(p, "doubleton", "second pair resistant")]
    assert (n[("BY=RM", "doubleton", "first pair resistant")], n[("BY=RM", "doubleton", "second pair resistant")], dt("BY=RM")) == (7, 5, 12)
    assert (n[("BY=YPS", "doubleton", "first pair resistant")], n[("BY=YPS", "doubleton", "second pair resistant")], dt("BY=YPS")) == (7, 3, 10)
    assert dt("BY=YJM") == 28 - 12 - 10 == 6


@pytest.mark.xfail(strict=True, reason=(
    "Text: with 50-kb windows, three loci (chrV near CAN1, chrX near ENT3/RSF2/VPS70, chrXIV near MKT1) "
    "were affected by more conditions than expected (8 or more distinct loci per bin). Binning File S10 "
    "positions (column position) in 50-kb bins from each chromosome start gives four bins with >= 8 loci: "
    "chrV 0-50 kb, chrX 650-700 kb, chrXII 650-700 kb and chrXV 600-650 kb; the chrXIV bin with MKT1 has fewer. "
    "The paper's bin boundaries are not given."))
def test_three_pleiotropic_50kb_bins(Q):
    b = collections.Counter((q["c"], q["p"] // 50000) for q in Q)
    assert sorted(k for k, v in b.items() if v >= 8) == [("chrV", 0), ("chrX", 13), ("chrXIV", 9)]


# ---------------------------------------------------------------- candidate genes -
def test_s11_is_the_135_singleton_loci(D, Q, CG):
    loci = {c["q1"] for c in CG}
    assert len(loci) == 135
    qid = {r: i for i, r in enumerate(D["qtl"]["id"])}
    by = collections.Counter(pname(D, c["sp"]) for c in {c["q1"]: c for c in CG}.values())
    assert by == {"BY": 42, "RM": 49, "YPS163": 26, "YJM789": 18}


def test_about_six_candidate_genes_per_locus(CG):
    # text: restricting to genes with nonsynonymous polymorphisms "reduced the average number to 6 per locus"
    n = collections.Counter(c["q1"] for c in CG)
    assert round(sum(n.values()) / len(n)) == 6


def test_candidate_genes_on_locus_chromosome_within_15kb(D, CG):
    G = D["genes"]
    for c in CG:
        if c["g"] < 0:
            continue
        q = D["qtl"]["p"][c["q1"]]
        assert G["c"][c["g"]] == D["qtl"]["c"][c["q1"]]
        assert max(0, G["a"][c["g"]] - q, q - G["e"][c["g"]]) <= 15000


def test_only_yer109c_missing_from_r64_5_1(CG):
    assert sorted({c["orf"] for c in CG if c["g"] < 0}) == ["YER109C"]


def test_tiers(D, ds):
    assert {t for t, d in zip(D["cand"]["tier"], D["cand"]["ds"]) if d == ds} == {"c"}


# ---------------------------------------------------------------- coordinates -----
def test_inside_r64_chromosomes_and_points(D, Q):
    for q in Q:
        assert 1 <= q["a"] <= q["p"] <= q["e"] <= D["chroms"][q["c"]]
        assert q["a"] == q["e"]


def test_known_loci_near_r64_genes(Q, D):
    """Loci the text places near CAN1, ENT3/RSF2/VPS70 and MKT1 sit at those genes' R64 positions."""
    G = rows(D["genes"]); pos = {g["n"]: g for g in G if g["n"]}
    def near(name, w=20000):
        g = pos[name]
        return sum(1 for q in Q if q["c"] == g["c"] and g["a"] - w <= q["p"] <= g["e"] + w)
    assert near("CAN1") >= 8 and near("VPS70") >= 8 and near("MKT1") >= 5


# ---------------------------------------------------------------- provenance ------
def test_provenance_resolves(D, ds):
    codes = {c for c, d in zip(D["sources"]["code"], D["sources"]["ds"]) if d == ds}
    assert codes == {"e12:S10", "e12:S11", "e12:S1dose", "e12:S1det"}
    for t in ("qtl", "causal", "crosses"):
        for pv, d in zip(D[t]["pv"], D[t]["ds"]):
            if d == ds:
                code, row = pv.split("#")
                assert code in codes and int(row) >= 2


def test_original_rows_kept(D):
    seen = collections.defaultdict(set)
    for ch in D["chroms"]:
        with open(os.path.join(ROOT, "docs", "v", f"e12_{ch}.json"), encoding="utf-8") as f:
            v = json.load(f)
        assert len(v["cols"]["e12:S10"]) == 17 and len(v["cols"]["e12:S11"]) == 7
        for code, rs in v["rows"].items():
            for rid, val in rs.items():
                assert len(val) == len(v["cols"][code])
                seen[code].add(tuple(val))
    assert len(seen["e12:S10"]) == 411 and len(seen["e12:S11"]) == 135
    with open(os.path.join(ROOT, "docs", "v", "e12_extra.json"), encoding="utf-8") as f:
        extra = json.load(f)
    assert len(extra["e12:S1dose"]["rows"]) == 13
