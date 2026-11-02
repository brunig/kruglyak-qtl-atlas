"""Validate the ingested Boocock et al. 2025 tables (single-cell eQTL, eLife 13:RP95566).

    python3 -I -m pytest -q tests/test_boocock2025.py      (after python3 -I ingest/boocock2025.py)

Numbers come from the article text (PMC11903034, quoted in sources/boocock2025/extracted/text_facts.tsv)
and the supplementary files. Disagreements would be strict xfails naming file and column.
"""
import collections
import csv
import io
import json
import os
import re
import tarfile
import zipfile

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "sources", "boocock2025")
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
    return D["datasets"]["id"].index("boocock2025")


@pytest.fixture(scope="module")
def L():
    with open(os.path.join(ROOT, "docs", "v", "b25_sceqtl.json"), encoding="utf-8") as f:
        return json.load(f)


@pytest.fixture(scope="module")
def Q(D, ds):
    return [r for r in rows(D["qtl"]) if r["ds"] == ds]


@pytest.fixture(scope="module")
def facts():
    with open(os.path.join(SRC, "extracted", "text_facts.tsv"), encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def lrows(L):
    return rows({k: v for k, v in L.items()})


def test_quotes_are_verbatim_in_the_article(facts):
    t = open(os.path.join(SRC, "article_pmc.xml"), encoding="utf-8").read()
    import html
    body = html.unescape(re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", t[t.find("<body"):t.find("<back")])))
    assert len(facts) >= 20
    for f in facts:
        assert f["quote"] in body, f["quote"]


def test_three_crosses_and_parents(D, ds):
    X = [r for r in rows(D["crosses"]) if r["ds"] == ds]
    P = D["parents"]["name"]
    assert {(x["code"], P[x["p1"]], P[x["p2"]]) for x in X} == {("A", "BY", "RM"), ("B", "YJM145", "YPS163"), ("C", "CBS2888", "YJM981")}


def test_local_eqtl_counts_per_cross(L):
    c = collections.Counter((x, l) for x, l in zip(L["xid"], L["loc"]))
    # text: 1031 local in cross A; 1193 in cross B and 721 in cross C
    assert (c[("b25-A", 1)], c[("b25-B", 1)], c[("b25-C", 1)]) == (1031, 1193, 721)


def test_distant_eqtl_counts_per_cross(L):
    c = collections.Counter((x, l) for x, l in zip(L["xid"], L["loc"]))
    # text: 1562 distant in cross A; 550 in cross B and 1126 in cross C (1626 in B and C together)
    assert (c[("b25-A", 0)], c[("b25-B", 0)], c[("b25-C", 0)]) == (1562, 550, 1126)


@pytest.mark.xfail(strict=True, reason=(
    "Text: '1626 distant eQTLs (550 in cross B and 1126 in cross C)'. The per-cross counts match Supplementary "
    "file 3 (sheets S2 and S3, rows with FDR (one-pot) < 0.05), but 550 + 1126 = 1676: the stated total is off by 50."))
def test_distant_total_in_new_crosses_1626(L):
    c = collections.Counter((x, l) for x, l in zip(L["xid"], L["loc"]))
    assert c[("b25-B", 0)] + c[("b25-C", 0)] == 1626


def test_local_total_in_new_crosses_1914(L):
    c = collections.Counter((x, l) for x, l in zip(L["xid"], L["loc"]))
    assert c[("b25-B", 1)] + c[("b25-C", 1)] == 1914


def test_distant_eqtl_are_on_another_chromosome(D, L):
    G = D["genes"]
    for g, c, l in zip(L["g"], L["c"], L["loc"]):
        if l == 0:
            assert G["c"][g] != c


def test_all_lazy_rows_are_significant(L):
    assert max(L["q"]) < 0.05


def test_hotspot_counts(Q, D):
    hs = [r for r in Q if r["k"].startswith("trans-eQTL hotspot")]
    c = collections.Counter(D["crosses"]["code"][r["x"]] for r in hs)
    # text: 12 hotspots in cross A; 13 in the new crosses (6 in B, 7 in C); 25 in total
    assert (c["A"], c["B"], c["C"], len(hs)) == (12, 6, 7, 25)


def test_hotspot_eqtl_counts_match_workbook_names(Q, L):
    hs = [r for r in Q if r["k"].startswith("trans-eQTL hotspot")]
    for r in hs:
        cx = r["id"].split("-")[2]
        n = sum(1 for x, l, c, h in zip(L["xid"], L["loc"], L["c"], L["hs"])
                if x == f"b25-{cx}" and l == 0 and h and h.split(":")[0] == r["c"]
                and r["a"] <= int(h.split(":")[1].split("-")[0]) and int(h.split("-")[1]) <= r["e"])
        assert n == r["cc"], r["id"]


def test_hotspots_are_the_supplementary_file_4_workbooks():
    with zipfile.ZipFile(os.path.join(SRC, "raw", "elife-95566-supp4-v1.zip")) as z:
        tgz = z.read("FileS4_VOR.tar.gz")
    with tarfile.open(fileobj=io.BytesIO(tgz), mode="r:gz") as t:
        names = [m for m in t.getnames() if re.fullmatch(r"fine_mapping/[ABC]/chr[IVX]+:\d+-\d+_\d+_finemapping\.xlsx", m)]
    assert len(names) == 25


def test_named_hotspot_genes(Q, D):
    G = D["genes"]
    got = {(r["id"].split("-")[2], r["c"], r["a"], r["e"]): (G["n"][r["g"]], r["gt"]) for r in Q if r["k"].startswith("trans") and r["g"] >= 0}
    assert got[("B", "chrVIII", 46887, 140660)] == ("GPA1", "x")
    assert got[("C", "chrX", 397734, 497167)] == ("CYR1", "c")
    assert got[("A", "chrX", 298300, 348016)] == ("ZAP1", "c")


def test_occupancy_qtl_20_unique(Q, D):
    occ = [r for r in Q if r["k"] == "cell-cycle occupancy QTL"]
    assert len(occ) == 27   # one row per QTL and stage in Supplementary file 7
    by = collections.defaultdict(list)
    for r in occ:
        by[D["crosses"]["code"][r["x"]]].append((r["c"], r["a"], r["e"]))
    merged = {}
    for cx, v in by.items():
        v.sort(); g = []
        for c, a, e in v:
            if g and g[-1][0] == c and a <= g[-1][2]:
                g[-1][2] = max(g[-1][2], e)
            else:
                g.append([c, a, e])
        merged[cx] = len(g)
    # text: "a total of 20 unique cell-cycle occupancy QTLs in the three crosses (4 for cross A, 10 for cross B, and 6 for cross C"
    assert (merged["A"], merged["B"], merged["C"]) == (4, 10, 6)


def test_occupancy_directions_match_text(Q, D):
    """GPA1: YJM145 allele more in G1; CYR1: CBS2888 allele more in G1 (positive = second parent)."""
    G, P = D["genes"], D["parents"]["name"]
    occ = {(G["n"][r["g"]], D["traits"]["cond"][r["t"]]): r for r in Q if r["k"] == "cell-cycle occupancy QTL" and r["g"] >= 0}
    assert P[occ[("GPA1", "G1")]["up"]] == "YJM145" and occ[("GPA1", "G1")]["ef"] < 0
    assert P[occ[("CYR1", "G1")]["up"]] == "CBS2888" and occ[("CYR1", "G1")]["ef"] < 0


def test_up_column_follows_effect_directions(L):
    second = {"b25-A": "RM", "b25-B": "YPS163", "b25-C": "YJM981"}
    first = {"b25-A": "BY", "b25-B": "YJM145", "b25-C": "CBS2888"}
    for x, ef, up in zip(L["xid"], L["ef"], L["up"]):
        assert up == (second[x] if ef > 0 else first[x])


def test_coordinates_inside_R64(Q, L):
    for r in Q:
        assert r["c"] in R64 and 0 <= r["a"] <= r["p"] <= r["e"] <= R64[r["c"]], r["id"]
    for c, a, e, p in zip(L["c"], L["a"], L["e"], L["p"]):
        assert c in R64 and 0 <= a <= p <= e <= R64[c]


def test_lazy_count_and_density(D, ds, L):
    lz = [r for r in rows(D["lazy"]) if r["ds"] == ds]
    assert len(lz) == 1 and lz[0]["name"] == "sceqtl" and lz[0]["n"] == len(L["id"])
    hb = D["hb"]
    re_c = collections.Counter(("local" if l else "distant", c, a // hb) for c, a, l in zip(L["c"], L["a"], L["loc"]))
    dn = {(r["k"], r["c"], r["b"]): r["n"] for r in rows(D["dens"]) if r["ds"] == ds}
    assert dn == dict(re_c)


def test_every_gene_resolves(L):
    assert min(L["g"]) >= 0


def test_provenance_resolves(D, ds, L):
    codes = {c for c, d in zip(D["sources"]["code"], D["sources"]["ds"]) if d == ds}
    for t in ("crosses", "qtl", "causal"):
        for r in rows(D[t]):
            if r.get("ds", ds) == ds and r["pv"]:
                assert r["pv"].split("#")[0] in codes
    assert all(p.split("#")[0] in codes for p in L["pv"])


def test_verbatim_rows_kept(L):
    seen = collections.Counter()
    for ch in R64:
        with open(os.path.join(ROOT, "docs", "v", f"b25_{ch}.json"), encoding="utf-8") as f:
            v = json.load(f)
        for code, rs in v["rows"].items():
            seen[code] += len(rs)
            for val in rs.values():
                assert len(val) == len(v["cols"][code])
    # every tested local transcript (supp2) and every distant eQTL (supp3) per cross
    assert seen["b25:S3A"] == 1562 and seen["b25:S3B"] == 550 and seen["b25:S3C"] == 1126
    assert seen["b25:S2A"] + seen["b25:S2B"] + seen["b25:S2C"] > 1031 + 1193 + 721
    assert seen["b25:S7A"] + seen["b25:S7B"] + seen["b25:S7C"] == 27
    assert seen["b25:S4"] == 25
    lazy_ids = set(L["id"])
    assert len(lazy_ids) == len(L["id"])
