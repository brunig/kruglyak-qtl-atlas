"""Validate the ingested Ehrenreich et al. 2010 tables against the paper's files.

    python3 -I -m pytest -q tests/test_ehrenreich2010.py      (after ./build.sh)

The readable sources are the Supplementary Information PDF and the article page's abstract
(the body is paywalled). Claims the files cannot reproduce are strict xfails naming the
file and column; claims nothing in the files can check are skipped with the reason.
"""
import collections
import html
import json
import os
import re

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "sources", "ehrenreich2010")
SI = os.path.join(SRC, "raw", "41586_2010_BFnature08923_MOESM285_ESM.pdf")
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
    return D["datasets"]["id"].index("ehrenreich2010")


@pytest.fixture(scope="module")
def Q(D, ds):
    return [r for r in rows(D["qtl"]) if r["ds"] == ds]


@pytest.fixture(scope="module")
def T(D, ds):
    return [r for r in rows(D["traits"]) if r["ds"] == ds]


@pytest.fixture(scope="module")
def CG(D, ds):
    return [r for r in rows(D["causal"]) if r["ds"] == ds]


@pytest.fixture(scope="module")
def si_text():
    pdfplumber = pytest.importorskip("pdfplumber")
    with pdfplumber.open(SI) as pdf:
        return re.sub(r"\s+", " ", " ".join((pg.extract_text() or "") for pg in pdf.pages))


@pytest.fixture(scope="module")
def facts():
    with open(os.path.join(SRC, "extracted", "text_facts.tsv")) as f:
        r = [l.rstrip("\n").split("\t") for l in f]
    return [dict(zip(r[0], x)) for x in r[1:]]


# ---------------------------------------------------------------- design ----------
def test_one_cross_BY_x_RM(D, ds):
    X = [r for r in rows(D["crosses"]) if r["ds"] == ds]
    assert len(X) == 1
    P = D["parents"]["name"]
    assert (P[X[0]["p1"]], P[X[0]["p2"]]) == ("BY", "RM")
    assert X[0]["n"] is None   # pools of more than 10^5 segregants, not counted


def test_17_chemicals_plus_mitochondrial_function(T):
    # abstract: "17 chemical resistance traits and mitochondrial function"; Table S1 lists 17 chemicals
    chem = [t for t in T if t["cat"] == "drug/chemical"]
    assert len(chem) == 17
    assert [t["name"] for t in T if t["cat"] != "drug/chemical"] == ["Mitochondrial function"]


def test_figure_s5_panels_match_table_s1(si_text):
    # Figures S5A-S5Q: one panel per chemical, 17 panels
    panels = re.findall(r"Figure S5([A-Q])\. Loci detected for ", si_text)
    assert panels == list("ABCDEFGHIJKLMNOPQ")


def test_condition_keys_link_to_bloom2019(D, ds):
    ck = {D["traits"]["ck"][t]: D["traits"]["ds"][t] for t in range(len(D["traits"]["id"]))}
    mine = {D["traits"]["ck"][t] for t in range(len(D["traits"]["id"])) if D["traits"]["ds"][t] == ds}
    if "bloom2019" in D["datasets"]["id"]:
        b = D["datasets"]["id"].index("bloom2019")
        bl = {D["traits"]["ck"][t] for t in range(len(D["traits"]["id"])) if D["traits"]["ds"][t] == b}
        assert {"neomycin 5mg/mL", "diamide 1.5mM"} <= (mine & bl)


# ---------------------------------------------------------------- loci ------------
def test_table_s2_two_known_4nqo_loci(Q, D):
    assert len(Q) == 2
    assert {D["traits"]["name"][q["t"]] for q in Q} == {"4-NQO 5 µg/mL"}
    assert sorted((q["c"], q["p"]) for q in Q) == [("chrXII", 203439), ("chrXIV", 470929)]
    assert all(q["m"] == "x" and q["lod"] is None and q["ef"] is None for q in Q)


@pytest.mark.xfail(strict=True, reason=(
    "Figure S4 legend: '14 4-NQO loci'. Only Table S2 of the Supplementary Information tabulates loci, and it has 2 rows "
    "(RAD5, MKT1); the other 12 are drawn in Figure S5A only and the main text is paywalled."))
def test_14_4nqo_loci(Q):
    assert len(Q) == 14


@pytest.mark.skip(reason="Abstract: 'some traits influenced by one major locus and others by at least 20 loci'. "
                         "Per-trait locus counts are only drawn (Figures S5A-S5Q); no table or readable text gives them.")
def test_loci_per_trait_range():
    pass


def test_table_s2_distance_is_peak_to_nearest_gene_end(CG, Q):
    with open(os.path.join(ROOT, "docs", "v", "e10_chrXII.json"), encoding="utf-8") as f:
        v12 = json.load(f)
    with open(os.path.join(ROOT, "docs", "v", "e10_chrXIV.json"), encoding="utf-8") as f:
        v14 = json.load(f)
    rowsv = list(v12["rows"]["e10:S2"].values()) + list(v14["rows"]["e10:S2"].values())
    for gene, chrom, span, peak, dist in rowsv:
        lo, hi = (int(v) for v in span.split("-"))
        assert dist == (lo - peak if peak < lo else peak - hi if peak > hi else 0), gene


def test_coordinates_inside_R64_and_ordered(Q, CG):
    for q in Q:
        assert 1 <= q["a"] <= q["p"] <= q["e"] <= R64[q["c"]]
    for c in CG:
        assert 1 <= c["a"] <= c["e"] <= R64[c["c"]]


def test_table_s2_gene_spans_inside_R64_genes(D, CG):
    """Build check: the source's RAD5 and MKT1 spans lie inside the R64-5-1 gene features of the same name."""
    G = D["genes"]
    for c in CG:
        if c["tier"] != "v":
            continue
        g = c["g"]
        assert G["n"][g] == c["name"] and G["c"][g] == c["c"]
        assert G["a"][g] <= c["a"] <= c["e"] <= G["e"][g], c["name"]


def test_candidate_genes_and_tiers(CG, Q, D):
    tiers = {(c["name"], c["tier"]) for c in CG}
    assert tiers == {("RAD5", "v"), ("MKT1", "v"), ("HAP1", "x")}
    ids = {r: i for i, r in enumerate(D["qtl"]["id"])}
    for q in Q:
        assert D["genes"]["n"][q["g"]] in ("RAD5", "MKT1") and q["gt"] == "v"
        assert any(c["q1"] == ids[q["id"]] for c in CG)


# ---------------------------------------------------------------- provenance ------
def test_text_facts_are_verbatim(facts, si_text):
    art = open(os.path.join(SRC, "article.html"), encoding="utf-8", errors="replace").read()
    art = re.sub(r"\s+", " ", html.unescape(re.sub(r"<[^>]+>", " ", art)))
    for f in facts:
        src = si_text if f["source"].endswith(".pdf") else art
        assert re.sub(r"\s+", " ", f["quote"]) in src, f["quote"][:60]


def test_provenance_resolves(D, ds):
    codes = {c for c, d in zip(D["sources"]["code"], D["sources"]["ds"]) if d == ds}
    assert all(c.startswith("e10:") for c in codes)
    assert D["sources"]["code"].count("e10:S2") == 1
    for t in ("crosses", "qtl", "causal"):
        for i, pv in enumerate(D[t]["pv"]):
            if t == "crosses" and D["crosses"]["ds"][i] != ds or t != "crosses" and D[t]["ds"][i] != ds:
                continue
            code, row = pv.split("#")
            assert code in codes and int(row) >= 2


def test_original_rows_kept(D):
    seen = collections.Counter()
    for ch in R64:
        with open(os.path.join(ROOT, "docs", "v", f"e10_{ch}.json"), encoding="utf-8") as f:
            v = json.load(f)
        for code, rs in v["rows"].items():
            for rid, val in rs.items():
                assert len(val) == len(v["cols"][code])
                seen[code] += 1
    assert seen == {"e10:S2": 4, "e10:TF": 1}   # 2 loci + their 2 causal rows; the HAP1 row
    with open(os.path.join(ROOT, "docs", "v", "e10_extra.json"), encoding="utf-8") as f:
        ex = json.load(f)
    assert len(ex["e10:S1"]["rows"]) == 17 and len(ex["e10:S2"]["rows"]) == 2
    assert {r[0] for r in ex["e10:SI"]["rows"]} == set(range(1, 23))   # every page's text kept
