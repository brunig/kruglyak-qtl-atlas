"""Validate the ingested Yvert et al. 2003 tables against the paper (abstract + supplements).

    python3 -I -m pytest -q tests/test_yvert2003.py     (after ingest)
"""
import collections
import json
import os
import re

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXT = os.path.join(ROOT, "sources", "yvert2003", "extracted")
R64 = {"chrI": 230218, "chrII": 813184, "chrIII": 316620, "chrIV": 1531933, "chrV": 576874, "chrVI": 270161,
       "chrVII": 1090940, "chrVIII": 562643, "chrIX": 439888, "chrX": 745751, "chrXI": 666816, "chrXII": 1078177,
       "chrXIII": 924431, "chrXIV": 784333, "chrXV": 1091291, "chrXVI": 948066}
# Supplementary Table 2 summary, as printed: group -> (messages, linkage bin, putative regulator)
SUMMARY = {1: (26, "II:390000", "Unknown"), 2: (272, "II:560000", "AMN1"), 3: (14, "II:710000", "Unknown"),
           4: (100, "III:100000", "LEU2"), 5: (36, "III:230000", "MAT"), 6: (19, "V:130000", "URA3"),
           7: (62, "VIII:130000", "GPA1"), 8: (114, "XII:680000", "HAP1"), 9: (34, "XII:1070000", "SIR3"),
           10: (21, "XIII:70000", "Unknown"), 11: (380, "XIV:503000", "Unknown"), 12: (174, "XV:180000", "Unknown"),
           13: (13, "XV:590000", "CAT5")}


def rows(tab):
    keys = list(tab)
    return [dict(zip(keys, v)) for v in zip(*(tab[k] for k in keys))]


def tsv(name):
    with open(os.path.join(EXT, name)) as f:
        lines = [l.rstrip("\n") for l in f]
    head = lines[1].split("\t")
    return [dict(zip(head, l.split("\t"))) for l in lines[2:] if l]


@pytest.fixture(scope="module")
def D():
    with open(os.path.join(ROOT, "docs", "data.json"), encoding="utf-8") as f:
        return json.load(f)


@pytest.fixture(scope="module")
def ds(D):
    return D["datasets"]["id"].index("yvert2003")


@pytest.fixture(scope="module")
def Q(D, ds):
    return [r for r in rows(D["qtl"]) if r["ds"] == ds]


@pytest.fixture(scope="module")
def L(D, ds):
    k = [i for i, d in enumerate(D["lazy"]["ds"]) if d == ds]
    assert len(k) == 1 and D["lazy"]["name"][k[0]] == "eqtl"
    with open(os.path.join(ROOT, "docs", D["lazy"]["file"][k[0]]), encoding="utf-8") as f:
        t = json.load(f)
    t["_n"] = D["lazy"]["n"][k[0]]
    return t


@pytest.fixture(scope="module")
def genes(D):
    return rows(D["genes"])


def test_quotes_are_verbatim():
    import html
    art = open(os.path.join(ROOT, "sources", "yvert2003", "article.html"), encoding="utf-8", errors="replace").read()
    abst = re.sub(r"\s+", " ", html.unescape(re.sub(r"<[^>]+>", " ", re.search(r'id="Abs1-content"[^>]*>(.*?)</div>', art, re.S).group(1))))
    pdfplumber = pytest.importorskip("pdfplumber")
    txt = ""
    for f in ("41588_2003_BFng1222_MOESM1_ESM.pdf", "41588_2003_BFng1222_MOESM2_ESM.pdf"):
        with pdfplumber.open(os.path.join(ROOT, "sources", "yvert2003", "raw", f)) as pdf:
            txt += " ".join((pg.extract_text() or "").replace("\n", " ") for pg in pdf.pages[:6]) + " "
    for r in tsv("text_facts.tsv"):
        assert r["quote"] in (abst if r["source"].startswith("article") else txt), r["quote"]


def test_citation_matches_page_metadata(D, ds):
    art = open(os.path.join(ROOT, "sources", "yvert2003", "article.html"), encoding="utf-8", errors="replace").read()
    assert '<meta name="citation_firstpage" content="57"' in art and '<meta name="citation_lastpage" content="64"' in art
    assert D["datasets"]["doi"][ds] == "10.1038/ng1222" and "35:57–64" in D["datasets"]["cite"][ds]


def test_one_cross_86_segregants(D, ds):
    X = [r for r in rows(D["crosses"]) if r["ds"] == ds]
    assert len(X) == 1 and X[0]["n"] == 86
    P = D["parents"]
    assert {P["name"][X[0]["p1"]], P["name"][X[0]["p2"]]} == {"BY", "RM"}


def test_13_groups_as_printed(Q):
    assert len(Q) == 13
    got = {int(re.search(r"group (\d+)", q["k"]).group(1)): (q["cc"], q["c"].replace("chr", "") + ":" + str(q["p"])) for q in Q}
    assert got == {g: (n, b) for g, (n, b, _) in SUMMARY.items()}


def test_22_bins_merged_into_13_groups():
    # "Twenty-two bins fit this criterion." / "Group 2 contains linkages from four bins, Group 11 from three bins,
    #  and Groups 3, 7, 8, and 12 from two bins each."
    bins = {g: 1 for g in SUMMARY}
    bins.update({2: 4, 11: 3, 3: 2, 7: 2, 8: 2, 12: 2})
    assert sum(bins.values()) == 22


def test_member_lists_match_summary_counts():
    c = collections.Counter(int(r["group"]) for r in tsv("group_members.tsv"))
    assert c == {g: n for g, (n, _, _) in SUMMARY.items()}
    assert sum(c.values()) == 1265


def test_group2_197_remaining_genes():
    # 12 (Figure 1 cluster) + 6 more DSE-related + 11 mitosis + 46 ribosome = 75; "the 197 Group 2 genes which are not DSE- or ribosome-related"
    assert SUMMARY[2][0] - (12 + 6 + 11 + 46) == 197


def test_group9_subtelomeric(genes):
    # "all but 2 of the 34 genes in Group 9 are located in the first or last 2% of chromosomes" (start sites on R64-5-1)
    gid = {}
    for i, g in enumerate(genes):
        gid.setdefault(g["id"], i)
    g9 = [r["orf"] for r in tsv("group_members.tsv") if r["group"] == "9"]
    assert len(g9) == 34
    inner = 0
    for o in g9:
        g = genes[gid[o]]
        start = g["a"] if g["s"] == "+" else g["e"]
        f = start / R64[g["c"]]
        inner += not (f <= 0.02 or f >= 0.98)
    assert inner == 2


def test_593_clusters():
    assert len(tsv("clusters.tsv")) == 593


def test_positions_inside_R64_points(Q, L):
    for q in Q:
        assert q["a"] == q["p"] == q["e"] and 0 < q["p"] <= R64[q["c"]]
    for c, a, p, e in zip(L["c"], L["a"], L["p"], L["e"]):
        assert a == p == e and 0 < p <= R64[c]


def test_named_regulators_near_their_bins(Q, genes):
    # build check: each named regulator on the group's chromosome within 50 kb of the printed bin label on R64-5-1
    byname = {}
    for i, g in enumerate(genes):
        if g["n"]:
            byname.setdefault(g["n"], g)
    pos = {g: (("chr" + b.split(":")[0]), int(b.split(":")[1])) for g, (_, b, _) in SUMMARY.items()}
    for grp, name in ((2, "AMN1"), (2, "MAK5"), (4, "LEU2"), (6, "URA3"), (7, "GPA1"), (8, "HAP1"), (9, "SIR3"), (13, "CAT5")):
        g = byname[name]; c, p = pos[grp]
        assert g["c"] == c and max(0, g["a"] - p, p - g["e"]) <= 50000, name


def test_candidate_tiers(D, ds, Q):
    C = [r for r in rows(D["causal"]) if r["ds"] == ds]
    tiers = {r["name"]: r["tier"] for r in C}
    assert tiers == {"AMN1": "x", "GPA1": "x", "MAK5": "c", "LEU2": "c", "URA3": "c", "HAP1": "c", "SIR3": "c", "CAT5": "c"}
    assert sum(1 for q in Q if q["g"] >= 0) == 7   # groups 1, 3, 5 (MAT), 10, 11, 12 name no single gene


def test_lazy_rows_and_density(D, ds, L):
    n = len(L["id"])
    assert n == L["_n"] == 1265
    assert all(v == 0 for v in L["loc"])          # cis-linkages were removed by the authors
    dens = collections.Counter((c, p // D["hb"]) for c, p in zip(L["c"], L["p"]))
    mine = {(c, b): k for c, b, k, d in zip(D["dens"]["c"], D["dens"]["b"], D["dens"]["n"], D["dens"]["ds"]) if d == ds}
    assert mine == dict(dens)


def test_unmatched_orfs_kept(L):
    assert sum(1 for g in L["g"] if g < 0) == 11


def test_no_invented_numbers(Q, L):
    assert all(q["lod"] is None and q["ef"] is None for q in Q)
    assert all(v is None for v in L["lod"]) and all(v is None for v in L["ef"])


def test_provenance_resolves(D, ds):
    codes = {c for c, d in zip(D["sources"]["code"], D["sources"]["ds"]) if d == ds}
    for t in ("qtl", "causal", "crosses"):
        for pv, d in zip(D[t]["pv"], D[t]["ds"] if "ds" in D[t] else [ds] * len(D[t]["pv"])):
            if d != ds:
                continue
            code, row = pv.split("#")
            assert code in codes and int(row) >= 3


def test_verbatim_rows_kept(D, ds):
    seen = collections.Counter()
    for ch in R64:
        with open(os.path.join(ROOT, "docs", "v", f"y03_{ch}.json"), encoding="utf-8") as f:
            v = json.load(f)
        for code, rs in v["rows"].items():
            seen[code] += len(rs)
    assert seen["y03:GM"] == 1265
    assert seen["y03:GS"] == 14     # 13 group rows + the printed continuation line under group 2
