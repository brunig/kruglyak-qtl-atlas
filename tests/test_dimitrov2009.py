"""Validate the ingested Dimitrov et al. 2009 tables against the paper (PMC full text) and themselves.

    python3 -I -m pytest -q tests/test_dimitrov2009.py      (after the ingest)

Numbers come from sources/dimitrov2009/raw/PMC2746160.html (text and Tables 1-4). Disagreements
inside the paper are strict xfails naming the table and the statement.
"""
import html
import json
import os
import re

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "sources", "dimitrov2009")
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
    return D["datasets"]["id"].index("dimitrov2009")


@pytest.fixture(scope="module")
def Q(D, ds):
    return [r for r in rows(D["qtl"]) if r["ds"] == ds]


@pytest.fixture(scope="module")
def CG(D, ds):
    return [r for r in rows(D["causal"]) if r["ds"] == ds]


@pytest.fixture(scope="module")
def TEXT():
    t = open(os.path.join(SRC, "raw", "PMC2746160.html"), encoding="utf-8", errors="replace").read()
    t = re.sub(r"<(script|style)[^>]*>.*?</\1>", " ", t, flags=re.S)
    return re.sub(r"\s+", " ", html.unescape(re.sub(r"<[^>]+>", " ", t)))


@pytest.fixture(scope="module")
def EXTRA():
    with open(os.path.join(ROOT, "docs", "v", "d09_extra.json"), encoding="utf-8") as f:
        return json.load(f)


def test_citation(D, ds, TEXT):
    assert D["datasets"]["doi"][ds] == "10.1534/genetics.109.104497"
    assert "Lazar N Dimitrov" in TEXT or 'content="Lazar N Dimitrov"' in open(os.path.join(SRC, "raw", "PMC2746160.html"), encoding="utf-8").read()


def test_quotes_verbatim(TEXT):
    with open(os.path.join(SRC, "extracted", "text_facts.tsv")) as f:
        facts = [l.rstrip("\n").split("\t") for l in f][1:]
    assert len(facts) == 28
    missing = [q for _, _, q in facts if q not in TEXT]
    assert not missing


def test_one_cross_122_segregants(D, ds):
    X = [r for r in rows(D["crosses"]) if r["ds"] == ds]
    assert len(X) == 1 and X[0]["n"] == 122
    P = D["parents"]
    assert {P["name"][X[0]["p1"]], P["name"][X[0]["p2"]]} == {"BY", "RM"}


def test_2957_markers_and_threshold(TEXT):
    assert "2957 polymorphic markers" in TEXT
    assert "linkage likelihood of 9.3 marks the 5% genomewide significance level" in TEXT


def test_three_qtl_on_xiv_and_xv(Q):
    assert len(Q) == 3
    assert [q["c"] for q in Q] == ["chrXIV", "chrXV", "chrXV"]


def test_windows_as_printed(Q):
    assert [(q["a"], q["e"]) for q in Q] == [(449639, 486861), (519776, 546197), (937036, 945781)]


def test_window_widths_match_text(Q, TEXT):
    # the legend gives 37.2, 26.4 and 8.7 kbp
    for q, w in zip(Q, ("37.2-kbp", "26.4-kbp", "8.7-kbp")):
        assert w in TEXT
        assert round((q["e"] - q["a"]) / 1000, 1) == float(w.split("-")[0])


def test_loci_2_and_3_far_apart(Q):
    # "likely to be unlinked genetically given the physical distance between them (~410 kbp)"
    assert 380000 < Q[2]["p"] - Q[1]["p"] < 440000


def test_coordinates_inside_r64(Q):
    for q in Q:
        assert 1 <= q["a"] <= q["p"] <= q["e"] <= R64[q["c"]]


def test_direction_by_allele_raises(D, Q):
    assert all(D["parents"]["name"][q["up"]] == "BY" for q in Q)


def test_four_genes_confirmed(CG, D):
    assert sorted(r["name"] for r in CG) == ["CAT5", "MIP1", "MKT1", "SAL1"]
    assert all(r["tier"] == "x" for r in CG)
    eff = {r["name"]: r["eff"] for r in CG}
    assert eff == {"SAL1": "increase", "CAT5": "increase", "MIP1": "increase", "MKT1": "decrease"}


def _gene(D, name):
    i = D["genes"]["n"].index(name)
    return D["genes"]["c"][i], D["genes"]["a"][i], D["genes"]["e"][i]


def test_genes_inside_their_windows_on_r64(D, Q):
    """Build check: the genes the paper places in loci 1 and 3 fall inside those windows on R64-5-1."""
    for name, q in (("MKT1", Q[0]), ("SAL1", Q[0]), ("MIP1", Q[2])):
        c, a, e = _gene(D, name)
        assert c == q["c"] and a >= q["a"] and e <= q["e"], name


def test_mip1_within_2kb_of_peak_window(D, Q):
    c, a, e = _gene(D, "MIP1")
    assert c == "chrXV" and a - 2000 <= Q[2]["e"] and e + 2000 >= Q[2]["a"]


def test_cat5_in_analysed_region_not_window(D, Q):
    """CAT5 lies right of the two-top-marker window but inside the 97-kbp region analysed (483,221-579,794)."""
    c, a, e = _gene(D, "CAT5")
    assert c == "chrXV" and a > Q[1]["e"] and 483221 <= a and e <= 579794


def test_table2_rates_from_counts(EXTRA):
    t = {r[0]: r[1:] for r in EXTRA["d09:T2"]["rows"]}
    two, one, rate = t["Two grandes"], t["Petite and grande"], t["Petite mutation rate (%)"]
    for g, p, r in zip(two, one, rate):
        assert round(100 * int(p) / (int(g) + int(p))) == int(r)


def test_table2_matches_text_rates(EXTRA, TEXT):
    t = {r[0]: r[1:] for r in EXTRA["d09:T2"]["rows"]}
    assert t["Petite mutation rate (%)"][:2] == ["5", "22"]
    assert "BY MKT1-30D strain was 5% while that of the BY MKT1-30G strain was 22%" in TEXT


@pytest.mark.xfail(strict=True, reason=(
    "The Results text gives the MKT1-30D vs MKT1-30G petite mutation rate difference as 'P = 0.01'; "
    "Table 2 (column UCC8241, row P-value) prints '<0.01a'."))
def test_mkt1_rate_pvalue_consistent(EXTRA):
    t = {r[0]: r[1:] for r in EXTRA["d09:T2"]["rows"]}
    assert t["P-value"][1] == "0.01"


def test_94_percent_explained(TEXT):
    assert round((79 - 5) / (79 - 0.6), 2) == 0.94
    assert "∼94%" in TEXT


def test_table3_rare_by_alleles(EXTRA):
    """SAL1 and CAT5: BY alleles only in S288C and BY4716; MIP1 also in W303 (text)."""
    head = EXTRA["d09:T3"]["rows"][0][1:]
    t = {r[0]: r[1:] for r in EXTRA["d09:T3"]["rows"][1:]}
    for g in ("SAL1", "CAT5"):
        assert [s for s, a in zip(head, t[g]) if a == "BY"] == ["S288C", "BY4716"]
    assert [s for s, a in zip(head, t["MIP1"]) if a == "BY"] == ["S288C", "BY4716", "W303"]


def test_tables_kept(EXTRA):
    assert [len(EXTRA[f"d09:T{k}"]["rows"]) for k in (1, 2, 3, 4)] == [34, 11, 5, 6]


def test_provenance(D, ds):
    codes = {c for c, d in zip(D["sources"]["code"], D["sources"]["ds"]) if d == ds}
    assert codes == {"d09:TF", "d09:T1", "d09:T2", "d09:T3", "d09:T4"}
    for t in ("crosses", "qtl", "causal"):
        for pv, d in zip(D[t]["pv"], D[t]["ds"]):
            if d == ds:
                code, row = pv.split("#")
                assert code in codes and int(row) >= 2


def test_verbatim_rows(Q, CG):
    seen = set()
    for ch in R64:
        with open(os.path.join(ROOT, "docs", "v", f"d09_{ch}.json"), encoding="utf-8") as f:
            v = json.load(f)
        for rid, val in v["rows"].get("d09:TF", {}).items():
            assert len(val) == 3
            seen.add(rid)
    assert {q["id"] for q in Q} | {c["id"] for c in CG} == seen
