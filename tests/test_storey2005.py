"""Validate the ingested Storey, Akey & Kruglyak 2005 tables against the paper and themselves.

    python3 -I -m pytest -q tests/test_storey2005.py      (after python3 -I ingest/storey2005.py)

Claims are asserted as the paper states them; where its own tables disagree the test is a strict
xfail naming the file and column.
"""
import collections
import json
import os
import re

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "sources", "storey2005")
R64 = {"chrI": 230218, "chrII": 813184, "chrIII": 316620, "chrIV": 1531933, "chrV": 576874, "chrVI": 270161,
       "chrVII": 1090940, "chrVIII": 562643, "chrIX": 439888, "chrX": 745751, "chrXI": 666816, "chrXII": 1078177,
       "chrXIII": 924431, "chrXIV": 784333, "chrXV": 1091291, "chrXVI": 948066}
ROMAN = list(R64)


def rows(tab):
    keys = list(tab)
    return [dict(zip(keys, v)) for v in zip(*(tab[k] for k in keys))]


def tsv(name):
    lines = open(os.path.join(SRC, "extracted", name)).read().splitlines()
    head = lines[1].split("\t")
    return [dict(zip(head, l.split("\t"))) for l in lines[2:]]


@pytest.fixture(scope="module")
def D():
    with open(os.path.join(ROOT, "docs", "data.json"), encoding="utf-8") as f:
        return json.load(f)


@pytest.fixture(scope="module")
def ds(D):
    return D["datasets"]["id"].index("storey2005")


@pytest.fixture(scope="module")
def E():
    with open(os.path.join(ROOT, "docs", "v", "s05_eqtl.json"), encoding="utf-8") as f:
        return json.load(f)


@pytest.fixture(scope="module")
def S1():
    return tsv("s1_table.tsv")


@pytest.fixture(scope="module")
def S2():
    return tsv("s2_table.tsv")


@pytest.fixture(scope="module")
def text():
    import html
    t = open(os.path.join(SRC, "article.html"), encoding="utf-8", errors="replace").read()
    t = re.sub(r"<(script|style)[^>]*>.*?</\1>", " ", t, flags=re.S)
    return " ".join(html.unescape(re.sub(r"<[^>]+>", " ", t)).split())


# ---------------------------------------------------------------- sources ----------
def test_quotes_verbatim(text):
    for q in tsv("text_facts.tsv"):
        assert " ".join(q["quote"].split()) in text, q["quote"][:60]


def test_extractions_name_current_pdfs():
    man = dict(l.split("\t")[0::2] for l in open(os.path.join(SRC, "MANIFEST.tsv")).read().splitlines()[1:])
    for name, pdf in (("s1_table.tsv", "journal.pbio.0030267.st001.pdf"), ("s2_table.tsv", "journal.pbio.0030267.st002.pdf")):
        first = open(os.path.join(SRC, "extracted", name)).readline()
        assert f"raw/{pdf} sha256={man[pdf]}" in first


# ---------------------------------------------------------------- design -----------
def test_one_cross_112_segregants(D, ds):
    X = [r for r in rows(D["crosses"]) if r["ds"] == ds]
    assert len(X) == 1 and X[0]["n"] == 112
    assert {D["parents"]["name"][X[0]["p1"]], D["parents"]["name"][X[0]["p2"]]} == {"BY", "RM"}


def test_170_traits_with_joint_linkage(S1):
    assert len(S1) == 170 and len({r["Trait"] for r in S1}) == 170


def test_58_cis_traits(S1):
    assert sum(1 for r in S1 if r["Cis"] != "0") == 58


@pytest.mark.xfail(strict=True, reason=(
    "The paper calls Table S1 the traits with joint linkage at FDR <= 10%, but 15 of its 170 rows print a "
    "joint q-value above 0.10 (largest 0.1049): journal.pbio.0030267.st001.pdf, column 'q-value Joint' "
    "(extracted/s1_table.tsv). Most defensible reading: rounding of a threshold applied before the printed values."))
def test_all_q_at_most_010(S1):
    assert max(float(r["q-value Joint"]) for r in S1) <= 0.10


def test_q_above_010_count_as_observed(S1):
    qs = [float(r["q-value Joint"]) for r in S1]
    assert sum(q > 0.10 for q in qs) == 15 and max(qs) == 0.1049


# ---------------------------------------------------------------- clusters (Table S2)
def test_ten_bins_with_three_or_more(S2):
    g = collections.Counter(r["Group"] for r in S2)
    assert len(g) == 10 and min(g.values()) >= 3


def test_named_cluster_sizes(S2):
    g = collections.Counter(r["Group"] for r in S2)
    fn = {r["Group"]: r["Putative Group Function"] for r in S2 if r["Putative Group Function"]}
    assert (g["1"], fn["1"]) == (12, "Mating")                                    # MAT (chr3) and GPA1 (chr8)
    assert (g["2"], fn["2"]) == (13, "Amino Acid and Mitochondria Metabolism")   # chr3 and chr13
    assert (g["3"], fn["3"]) == (10, "Mitochondrial TCA")                         # opposite ends of chr15
    assert (g["5"], fn["5"]) == (5, "Response To Stress")                         # chr6 and chr10


def _chr15_groups(S2):
    return {r["Group"] for r in S2 if ("15", "172654") in ((r["Locus 1 Chr"], r["Locus 1 Pos"]), (r["Locus 2 Chr"], r["Locus 2 Pos"]))}


@pytest.mark.xfail(strict=True, reason=(
    "Text: of the 10 bins 'there are three that share a common linkage to the exact same position on Chromosome 15'. "
    "Table S2 (journal.pbio.0030267.st002.pdf, columns 'Locus 1/2 Chr' and 'Locus 1/2 Pos') has four groups with a "
    "locus at chr15 172654 (groups 3, 4, 6 and 7 by the extraction)."))
def test_three_bins_share_chr15_172654(S2):
    assert len(_chr15_groups(S2)) == 3


def test_chr15_groups_as_observed(S2):
    assert len(_chr15_groups(S2)) == 4


def test_mating_cluster_on_chr3_and_chr8(S2):
    for r in S2:
        if r["Group"] == "1":
            assert {r["Locus 1 Chr"], r["Locus 2 Chr"]} == {"3", "8"}


def _s2_vs_s1(S1, S2):
    """Table S2 rows that disagree with Table S1 on the locus pair, the q-value or whether a locus is cis
    (S2's Cis column is yes/no; S1's says which locus is cis)."""
    s1 = {r["Trait"]: r for r in S1}
    bad = []
    for r in S2:
        a = s1.get(r["Trait"])
        pair = lambda x: {(x["Locus 1 Chr"], x["Locus 1 Pos"]), (x["Locus 2 Chr"], x["Locus 2 Pos"])}
        if a is None or pair(a) != pair(r) or a["q-value Joint"] != r["Joint FDR"] or (a["Cis"] != "0") != (r["Cis"] != "0"):
            bad.append(r["Trait"])
    return bad


@pytest.mark.xfail(strict=True, reason=(
    "Table S2's row for YFL059W (journal.pbio.0030267.st002.pdf) prints its chr6 locus at position 52 with Cis 0; "
    "Table S1 (st001.pdf) prints the same trait's chr6 locus at 5852 with Cis 1 (cis to locus 1). Every other Table S2 "
    "row agrees with Table S1."))
def test_every_table_s2_row_matches_table_s1(S1, S2):
    assert not _s2_vs_s1(S1, S2)


def test_table_s2_mismatch_is_only_yfl059w(S1, S2):
    assert _s2_vs_s1(S1, S2) == ["YFL059W"]


# ---------------------------------------------------------------- lazy table ---------
def test_lazy_count(D, ds, E):
    L = D["lazy"]
    k = [i for i in range(len(L["ds"])) if L["ds"][i] == ds]
    assert len(k) == 1 and L["name"][k[0]] == "eqtl" and L["n"][k[0]] == len(E["id"]) == 340


def test_two_loci_per_trait(E):
    c = collections.Counter(pv for pv in E["pv"])
    assert set(c.values()) == {2} and len(c) == 170
    assert collections.Counter(E["lk"]) == {1: 170, 2: 170}


def test_points_inside_r64(E):
    for c, a, e, p in zip(E["c"], E["a"], E["e"], E["p"]):
        assert a == e == p and 1 <= p <= R64[c]


def test_local_rows_follow_cis_flag(E, S1):
    assert sum(E["loc"]) == 58


def test_cis_locus_on_trait_chromosome(S1):
    for r in S1:
        if r["Cis"] != "0":
            assert r[f"Locus {r['Cis']} Chr"] == r["Chr"], r["Trait"]


def test_cis_locus_near_trait_gene(S1):
    far = []
    for r in S1:
        if r["Cis"] != "0":
            p = int(r[f"Locus {r['Cis']} Pos"]); a, b = sorted((int(r["Start"]), int(r["Stop"])))
            if max(0, a - p, p - b) > 50000:
                far.append(r["Trait"])
    assert not far


def test_genes_resolve_except_retired_orfs(D, E):
    gids = {D["genes"]["id"][g] for g in E["g"] if g >= 0}
    assert len(gids) == 163
    assert sum(1 for g in E["g"] if g < 0) == 14   # 7 trait ORFs absent from R64-5-1, two rows each


def test_dens_recount(D, ds, E):
    hb = D["hb"]
    want = collections.Counter(("local" if l else "distant", c, p // hb) for l, c, p in zip(E["loc"], E["c"], E["p"]))
    got = collections.Counter({(k, c, b): n for d, k, c, b, n in zip(D["dens"]["ds"], D["dens"]["k"], D["dens"]["c"], D["dens"]["b"], D["dens"]["n"]) if d == ds})
    assert got == want


def test_no_qtl_rows_no_invented_effects(D, ds, E):
    assert sum(1 for d in D["qtl"]["ds"] if d == ds) == 0
    assert all(v is None for v in E["lod"]) and all(v is None for v in E["ef"])


# ---------------------------------------------------------------- provenance ---------
def test_provenance_resolves(D, ds, E):
    codes = {c for c, d in zip(D["sources"]["code"], D["sources"]["ds"]) if d == ds}
    n = len(tsv("s1_table.tsv"))
    for pv in E["pv"]:
        code, row = pv.split("#")
        assert code in codes and 2 <= int(row) <= n + 1
    for pv, d in zip(D["crosses"]["pv"], D["crosses"]["ds"]):
        if d == ds:
            code, row = pv.split("#")
            assert code in codes and int(row) >= 2


def test_verbatim_rows_kept(E):
    seen = 0
    for c in R64:
        with open(os.path.join(ROOT, "docs", "v", f"s05_{c}.json"), encoding="utf-8") as f:
            v = json.load(f)
        for rid, row in v["rows"]["s05:S1"].items():
            assert len(row) == len(v["cols"]["s05:S1"])
            seen += 1
    assert seen == 340
    with open(os.path.join(ROOT, "docs", "v", "s05_extra.json"), encoding="utf-8") as f:
        x = json.load(f)
    assert len(x["s05:S1"]["rows"]) == 170 and len(x["s05:S2"]["rows"]) == 61
