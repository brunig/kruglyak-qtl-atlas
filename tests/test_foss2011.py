"""Validate the Foss et al. 2011 ingest against the paper's own numbers and against itself.

    python3 -I -m pytest -q tests/test_foss2011.py      (after ./build.sh or ingest/foss2011.py)
"""
import collections
import json
import os
import re
import html

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "sources", "foss2011")
RAW = os.path.join(SRC, "raw")
R64 = {"chrI": 230218, "chrII": 813184, "chrIII": 316620, "chrIV": 1531933, "chrV": 576874, "chrVI": 270161,
       "chrVII": 1090940, "chrVIII": 562643, "chrIX": 439888, "chrX": 745751, "chrXI": 666816, "chrXII": 1078177,
       "chrXIII": 924431, "chrXIV": 784333, "chrXV": 1091291, "chrXVI": 948066}


@pytest.fixture(scope="module")
def D():
    with open(os.path.join(ROOT, "docs", "data.json"), encoding="utf-8") as f:
        return json.load(f)


@pytest.fixture(scope="module")
def ds(D):
    return D["datasets"]["id"].index("foss2011")


def lazy(D, ds, name):
    L = D["lazy"]
    k = [i for i in range(len(L["ds"])) if L["ds"][i] == ds and L["name"][i] == name][0]
    with open(os.path.join(ROOT, "docs", L["file"][k]), encoding="utf-8") as f:
        return L["n"][k], json.load(f)


@pytest.fixture(scope="module")
def P(D, ds):
    return lazy(D, ds, "mspqtl")


@pytest.fixture(scope="module")
def E(D, ds):
    return lazy(D, ds, "eqtl")


def s4():
    import xlrd
    sh = xlrd.open_workbook(os.path.join(RAW, "journal.pbio.1001144.s009.xls")).sheet_by_index(0)
    return [sh.row_values(i) for i in range(1, sh.nrows)]


def article_text():
    t = open(os.path.join(SRC, "article.html"), encoding="utf-8", errors="replace").read()
    t = html.unescape(re.sub(r"<[^>]+>", " ", t))
    return re.sub(r"\s+", " ", t)


# ---------------------------------------------------------------- the paper's numbers --------
def test_one_cross_95_segregants(D, ds):
    X = [i for i, d in enumerate(D["crosses"]["ds"]) if d == ds]
    assert len(X) == 1 and D["crosses"]["n"][X[0]] == 95
    assert sorted(D["parents"]["name"][i] for i, d in enumerate(D["parents"]["ds"]) if d == ds) == ["BY", "RM"]


def test_linkages_179_proteins_342_transcripts(P, E):
    assert P[0] == len(P[1]["id"]) == 179
    assert E[0] == len(E[1]["id"]) == 342


def test_125_proteins_200_transcripts_with_a_linkage(P, E):
    assert len(set(P[1]["g"])) == 125
    assert len(set(E[1]["g"])) == 200


@pytest.mark.xfail(strict=True, reason=(
    "Text: 'At a 5% FDR, we mapped 49 and 97 loci' (proteins, transcripts). Table S4 (s009.xls, column marker) has 73 "
    "distinct markers for protein linkages and 156 for transcript linkages; the rule that merged markers into loci is not "
    "stated, so the locus counts cannot be reproduced."))
def test_49_and_97_loci():
    rows = s4()
    assert (len({r[4] for r in rows if r[0] == "peptide"}), len({r[4] for r in rows if r[0] == "transcript"})) == (49, 97)


def test_distinct_markers_as_observed():
    rows = s4()
    assert (len({r[4] for r in rows if r[0] == "peptide"}), len({r[4] for r in rows if r[0] == "transcript"})) == (73, 156)


def test_354_proteins_95_segregants_in_table_s1():
    import xlrd
    sh = xlrd.open_workbook(os.path.join(RAW, "journal.pbio.1001144.s006.xls")).sheet_by_index(0)
    assert sh.nrows - 1 == 354 and sh.ncols - 3 == 95


def test_94_of_354_correlations_significant_and_cutoffs():
    import statistics
    import xlrd
    sh = xlrd.open_workbook(os.path.join(RAW, "journal.pbio.1001144.s008.xls")).sheet_by_index(0)
    rows = [sh.row_values(i) for i in range(1, sh.nrows)]
    assert len(rows) == 354 and sum(1 for r in rows if r[5] == 1) == 94
    assert round(statistics.mean(r[3] for r in rows), 2) == 0.24 and round(statistics.mean(r[4] for r in rows), 2) == 0.31


def test_local_counts_consistent_with_paper(P, E):
    # 1% FDR (a subset of these 5% FDR linkages): 22 transcripts and 3 proteins act in cis
    assert len({E[1]["g"][i] for i in range(len(E[1]["g"])) if E[1]["loc"][i] == 1}) >= 22
    assert len({P[1]["g"][i] for i in range(len(P[1]["g"])) if P[1]["loc"][i] == 1}) >= 3


def test_hotspots_at_leu2_and_amn1(D):
    rows = s4()
    c = collections.Counter((int(r[5]), int(r[6])) for r in rows if r[0] == "peptide")
    G = D["genes"]
    def gene(n):
        i = G["n"].index(n)
        return G["c"][i], G["a"][i], G["e"][i]
    for name, chrom_no in (("LEU2", 3), ("AMN1", 2)):
        gc, ga, ge = gene(name)
        near = sum(n for (cc, p), n in c.items() if cc == chrom_no and ga - 20000 <= p <= ge + 20000)
        assert near >= 8, name


# ---------------------------------------------------------------- self-consistency ------------
def test_positions_inside_r64_points(P, E):
    for L in (P[1], E[1]):
        for i in range(len(L["id"])):
            assert L["c"][i] in R64 and 1 <= L["p"][i] <= R64[L["c"][i]]
            assert L["a"][i] == L["p"][i] == L["e"][i]


def test_marker_names_match_positions():
    for r in s4():
        parts = r[4].split("_")
        assert int(parts[-2]) == int(r[5]) and int(parts[-1]) == int(r[6])


def test_every_gene_resolves(P, E):
    assert all(g >= 0 for g in P[1]["g"] + E[1]["g"])


def test_no_invented_effects(P, E):
    assert all(v is None for v in P[1]["ef"] + E[1]["ef"] + P[1]["lod"] + E[1]["lod"])


def test_density_recount(D, ds, P, E):
    got = collections.Counter()
    for label, L in (("protein", P[1]), ("mRNA", E[1])):
        for i in range(len(L["id"])):
            got[(label + (", local" if L["loc"][i] else ", distant"), L["c"][i], L["p"][i] // D["hb"])] += 1
    Dn = D["dens"]
    have = collections.Counter({(Dn["k"][i], Dn["c"][i], Dn["b"][i]): Dn["n"][i] for i in range(len(Dn["ds"])) if Dn["ds"][i] == ds})
    assert got == have


def test_quotes_verbatim():
    t = article_text()
    with open(os.path.join(SRC, "extracted", "text_facts.tsv")) as f:
        lines = [l.rstrip("\n").split("\t") for l in f if not l.startswith("#")][1:]
    for k, q in lines:
        assert re.sub(r"\s+", " ", q) in t, k


def test_provenance_resolves(D, ds, P, E):
    codes = {D["sources"]["code"][i] for i in range(len(D["sources"]["code"])) if D["sources"]["ds"][i] == ds}
    assert codes == {"f11:S1", "f11:S3", "f11:S4", "f11:TF"}
    for L in (P[1], E[1]):
        for pv in L["pv"]:
            code, row = pv.split("#")
            assert code == "f11:S4" and 2 <= int(row) <= 522


def test_verbatim_rows_kept(P, E):
    ids = set(P[1]["id"]) | set(E[1]["id"])
    seen = set()
    for c in R64:
        with open(os.path.join(ROOT, "docs", "v", f"f11_{c}.json"), encoding="utf-8") as f:
            v = json.load(f)
        for rid, row in v["rows"]["f11:S4"].items():
            assert len(row) == len(v["cols"]["f11:S4"]) == 7
            seen.add(rid)
    assert seen == ids
    with open(os.path.join(ROOT, "docs", "v", "f11_extra.json"), encoding="utf-8") as f:
        x = json.load(f)
    assert len(x["f11:S1"]["rows"]) == 354 and len(x["f11:S3"]["rows"]) == 354
