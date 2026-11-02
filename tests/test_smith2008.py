"""Validate the ingested Smith & Kruglyak 2008 tables against the paper and against themselves.

    python3 -I -m pytest -q tests/test_smith2008.py      (after python3 -I ingest/smith2008.py)

Numbers quoted from the article are checked against sources/smith2008/extracted/text_facts.tsv (each
quote verbatim in the page) and against the ingested rows. Discrepancies would be strict xfails naming
file and column.
"""
import collections
import html
import json
import os
import re

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "sources", "smith2008")
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
    return D["datasets"]["id"].index("smith2008")


@pytest.fixture(scope="module")
def L():
    with open(os.path.join(ROOT, "docs", "v", "s08_eqtl.json"), encoding="utf-8") as f:
        return json.load(f)


@pytest.fixture(scope="module")
def clean(L):
    return [r for r in rows(L) if r["st"] is None]


@pytest.fixture(scope="module")
def Q(D, ds):
    return [r for r in rows(D["qtl"]) if r["ds"] == ds]


@pytest.fixture(scope="module")
def facts():
    with open(os.path.join(SRC, "extracted", "text_facts.tsv"), encoding="utf-8") as f:
        lines = [l.rstrip("\n").split("\t") for l in f][2:]
    return {r[0]: r[2] for r in lines}


def page_text():
    t = open(os.path.join(SRC, "article.html"), encoding="utf-8", errors="replace").read()
    t = re.sub(r"<(script|style)[^>]*>.*?</\1>", " ", t, flags=re.S)
    return re.sub(r"\s+", " ", html.unescape(re.sub(r"<[^>]+>", " ", t)))


def test_quotes_verbatim(facts):
    text = page_text()
    assert len(facts) == 29
    for k, q in facts.items():
        assert q in text, k


def test_design(D, ds, facts):
    X = [r for r in rows(D["crosses"]) if r["ds"] == ds]
    assert len(X) == 1 and X[0]["n"] == 109
    P = D["parents"]
    assert {P["name"][X[0]["p1"]], P["name"][X[0]["p2"]]} == {"BY", "RM"}
    assert "109 segregant strains" in facts["design"]


def test_linkage_counts_match_text(clean, facts):
    """3,997 glucose, 3,489 ethanol and 1,555 gxe linkages: Table S4 rows without probe polymorphisms."""
    c = collections.Counter(r["env"] for r in clean)
    assert (c["glucose"], c["ethanol"], c["gxe"]) == (3997, 3489, 1555)
    assert "3,997 and 3,489" in facts["linkages_within"] and "total of 1,555 linkages" in facts["gxe_total"]


def test_gxe_transcripts_1382(D, L, clean):
    ids = {r["pv"] for r in clean if r["env"] == "gxe"}
    assert len(ids) == 1555
    # count distinct transcripts through the verbatim rows (ORF column), since some ORFs are not in R64-5-1
    orfs = set()
    for ch in R64:
        v = json.load(open(os.path.join(ROOT, "docs", "v", f"s08_{ch}.json"), encoding="utf-8"))
        for lid, row in v["rows"]["s08:S4"].items():
            if row[2] == "gxeQTL" and row[1] == "no":
                orfs.add(row[0])
    assert len(orfs) == 1382


def test_probe_polymorphism_rows_flagged(L):
    st = collections.Counter(L["st"])
    assert st[None] == 9041 and len(L["id"]) == 11392
    assert st["probe polymorphism: excluded by the paper"] == 2271


def test_local_counts(clean, facts):
    loc = collections.Counter((r["env"], r["loc"]) for r in clean)
    assert loc[("gxe", 1)] == 172 and loc[("gxe", 0)] == 1383
    assert "172/1,555 (11%)" in facts["local_gxe"]
    # "22%-25% of all linkages were local" within conditions
    for env in ("glucose", "ethanol"):
        frac = loc[(env, 1)] / (loc[(env, 1)] + loc[(env, 0)])
        assert 0.215 <= frac <= 0.255, (env, frac)


def test_lod_thresholds(clean, facts):
    m = {e: min(r["lod"] for r in clean if r["env"] == e) for e in ("glucose", "ethanol", "gxe")}
    assert m["glucose"] >= 3.25 and m["ethanol"] >= 3.3 and m["gxe"] >= 3.7
    assert "3.25, 3.3, and 3.7" in facts["thresholds"]


def test_peaks_13_13_15(Q, facts):
    c = collections.Counter(r["k"].split()[2][:3] for r in Q)
    assert (c["glu"], c["eth"], c["gxe"]) == (13, 13, 15)
    assert "13, 13, and 15 peaks" in facts["peaks"]


def test_peak_linkage_counts_match_table_s1(Q):
    import xlrd
    s = xlrd.open_workbook(os.path.join(SRC, "raw", "journal.pbio.0060083.st001.xls")).sheet_by_index(0)
    tot = {s.cell_value(i, 0): int(s.cell_value(i, 3)) for i in range(1, s.nrows) if isinstance(s.cell_value(i, 1), float)}
    assert len(tot) == 41
    for r in Q:
        assert r["cc"] == tot[r["k"].split()[2]], r["k"]


def test_named_peak_counts(Q, facts):
    cc = {r["k"].split()[2]: r["cc"] for r in Q}
    assert cc["gxe12"] == 372 and "372 transcripts" in facts["ira2_peak"]
    assert cc["eth11"] in (386, 387)   # text: 386 transcripts; Table S1 total linkages 387
    assert cc["glu12"] == 1159 and cc["eth10"] == 410 and "(1,159 and 410" in facts["ira2_within"]
    assert cc["glu7"] == 49 and "49 transcripts" in facts["gpa1"]


def _distant_gxe_in_peaks():
    v = {}
    for ch in R64:
        v.update(json.load(open(os.path.join(ROOT, "docs", "v", f"s08_{ch}.json"), encoding="utf-8"))["rows"]["s08:S4"])
    return sum(1 for lid, row in v.items() if row[2] == "gxeQTL" and row[1] == "no" and row[10] == "distant" and str(row[9]).startswith("gxe"))


@pytest.mark.xfail(strict=True, reason=(
    "Text: 81% (1,126/1,383) of distant gxe QTL fall in significant bins. Table S4 (st004.xls), column "
    "'distant peak', assigns 1,123 polymorphism-free distant gxe linkages to the 15 gxe peaks, and Table S1 "
    "(st001.xls), column 'Total linkages', sums to the same 1,123 over gxe1-gxe15."))
def test_distant_gxe_in_peaks_1126(facts):
    assert "1,126/1,383" in facts["distant_bins"]
    assert _distant_gxe_in_peaks() == 1126


def test_distant_gxe_in_peaks_observed():
    assert _distant_gxe_in_peaks() == 1123


def test_coordinates(L, Q):
    beyond = []
    for c, a, p, e in zip(L["c"], L["a"], L["p"], L["e"]):
        assert c in R64 and 0 <= a <= p <= e   # a support interval can start at the chromosome start
        assert p <= R64[c]
        if e > R64[c]:
            beyond.append((c, e))
    # one right-hand marker, chrVII:1,090,946, lies 6 bp past the R64 chrVII end (1,090,940): kept as printed
    assert set(beyond) == {("chrVII", 1090946)} and len(beyond) == 8
    for r in Q:
        assert r["c"] in R64 and r["a"] <= r["p"] <= r["e"] <= R64[r["c"]]


def test_effect_sign_and_up(L):
    n = 0
    for ef, up in zip(L["ef"], L["up"]):
        if ef is None:
            continue
        assert up == ("RM" if ef > 0 else "BY" if ef < 0 else None)
        n += 1
    assert n > 11000


def test_candidates(D, ds, Q):
    C = [r for r in rows(D["causal"]) if r["ds"] == ds]
    tiers = {r["name"]: r["tier"] for r in C if r["tier"] == "x"}
    assert {"AMN1", "GPA1", "HAP1", "IRA2", "DIG1", "MKT1", "MAT"} >= set(tiers) and {"IRA2", "DIG1", "MKT1"} <= set(tiers)
    G = D["genes"]
    qid = {r["id"]: r for r in Q}
    ids = D["qtl"]["id"]
    for r in C:
        if r["g"] < 0:
            continue
        q = qid[ids[r["q1"]]]
        assert G["c"][r["g"]] == q["c"]
        # every named gene lies within 25 kb of the span of its peak's linkages
        d = max(0, q["a"] - G["e"][r["g"]], G["a"][r["g"]] - q["e"])
        assert d <= 25000, (r["name"], d)


def test_lazy_and_density(D, ds, L, clean):
    lz = [r for r in rows(D["lazy"]) if r["ds"] == ds]
    assert len(lz) == 1 and lz[0]["n"] == len(L["id"]) == 11392
    dens = collections.Counter()
    for r in clean:
        for b in range(r["a"] // 50000, r["e"] // 50000 + 1):
            dens[(r["env"], r["c"], b)] += 1
    got = {(k, c, b): n for k, c, b, n, d in zip(D["dens"]["k"], D["dens"]["c"], D["dens"]["b"], D["dens"]["n"], D["dens"]["ds"]) if d == ds}
    assert got == dict(dens)


def test_genes_resolve(L):
    assert sum(1 for g in L["g"] if g >= 0) / len(L["g"]) > 0.99


def test_provenance_and_verbatim(D, ds, L, Q):
    codes = {c for c, d in zip(D["sources"]["code"], D["sources"]["ds"]) if d == ds}
    for pv in L["pv"] + [r["pv"] for r in Q]:
        code, row = pv.split("#")
        assert code in codes and int(row) >= 2
    seen = collections.Counter()
    for ch in R64:
        v = json.load(open(os.path.join(ROOT, "docs", "v", f"s08_{ch}.json"), encoding="utf-8"))
        for code, rs in v["rows"].items():
            seen[code] += len(rs)
            assert all(len(x) == len(v["cols"][code]) for x in rs.values())
    assert seen == {"s08:S4": 11392, "s08:S1": 41}
