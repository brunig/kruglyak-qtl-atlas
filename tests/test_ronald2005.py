"""Validate the ingested Ronald et al. 2005 tables against the paper's text and its own files.

    python3 -I -m pytest -q tests/test_ronald2005.py      (after the ingest)

Paper claims are asserted as stated; where the files disagree, the test is a strict xfail naming
file and column. The article text is the tag-stripped page in sources/ronald2005/article.html.
"""
import html
import json
import math
import os
import re

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "sources", "ronald2005")
RAW = os.path.join(SRC, "raw")
P = "journal.pgen.0010025."
R64 = {"chrI": 230218, "chrII": 813184, "chrIII": 316620, "chrIV": 1531933, "chrV": 576874, "chrVI": 270161,
       "chrVII": 1090940, "chrVIII": 562643, "chrIX": 439888, "chrX": 745751, "chrXI": 666816, "chrXII": 1078177,
       "chrXIII": 924431, "chrXIV": 784333, "chrXV": 1091291, "chrXVI": 948066}


def rows_of(name):
    with open(os.path.join(RAW, P + name), encoding="utf-8") as f:
        lines = [l.rstrip("\r\n") for l in f]
    return lines[0].split("\t"), [l.split("\t") for l in lines[1:] if l.strip()]


@pytest.fixture(scope="module")
def D():
    with open(os.path.join(ROOT, "docs", "data.json"), encoding="utf-8") as f:
        return json.load(f)


@pytest.fixture(scope="module")
def ds(D):
    return D["datasets"]["id"].index("ronald2005")


@pytest.fixture(scope="module")
def L(D, ds):
    k = [i for i, d in enumerate(D["lazy"]["ds"]) if d == ds][0]
    with open(os.path.join(ROOT, "docs", D["lazy"]["file"][k]), encoding="utf-8") as f:
        t = json.load(f)
    t["_n"] = D["lazy"]["n"][k]
    t["_name"] = D["lazy"]["name"][k]
    return t


@pytest.fixture(scope="module")
def S1():
    return rows_of("st001.txt")


@pytest.fixture(scope="module")
def S2():
    return rows_of("st002.txt")


@pytest.fixture(scope="module")
def text():
    t = open(os.path.join(SRC, "article.html"), encoding="utf-8", errors="replace").read()
    m = re.search(r'<div class="article-text"[^>]*>(.*)', t, re.S)
    t = m.group(1) if m else t
    t = re.sub(r"<(script|style)[^>]*>.*?</\1>", " ", t, flags=re.S)
    return re.sub(r"\s+", " ", html.unescape(re.sub(r"<[^>]+>", " ", t)))


def test_quotes_verbatim(text):
    with open(os.path.join(SRC, "extracted", "text_facts.tsv")) as f:
        lines = [l.rstrip("\n").split("\t") for l in f if not l.startswith("#")][1:]
    assert len(lines) == 19
    for k, sec, q in lines:
        assert q in text, k


def test_design(D, ds):
    X = [i for i, d in enumerate(D["crosses"]["ds"]) if d == ds]
    assert len(X) == 1 and D["crosses"]["n"][X[0]] == 112
    names = {D["parents"]["name"][p] for p, d in enumerate(D["parents"]["ds"]) if d == ds}
    assert names == {"BY", "RM"}
    assert D["datasets"]["fdr"][ds] == 0.05


def test_5727_tested_1428_self_linkages(S1):
    head, rows = S1
    assert len(rows) == 5727
    assert sum(1 for r in rows if r[2] == "1") == 1428


def test_self_linkage_threshold_p_0012(S1):
    # every self-linkage has p < 0.012 (minus_ln_pvalue > -ln 0.012)
    sl = [float(r[6]) for r in S1[1] if r[2] == "1"]
    assert min(sl) >= -math.log(0.0122)
    assert all(r[3] == "NA" for r in S1[1] if r[2] == "0")


def test_446_with_effect_above_1_2(S1):
    assert sum(1 for r in S1[1] if r[2] == "1" and abs(float(r[3])) > math.log2(1.2)) == 446


def test_5182_aligned_1233_self(S1):
    rows = S1[1]
    assert sum(1 for r in rows if r[1] == "1") == 5182
    assert sum(1 for r in rows if r[1] == "1" and r[2] == "1") == 1233


def test_ase_77_self_16_controls(S2):
    rows = S2[1]
    sl = [r for r in rows if r[1] == "self.link"]
    ctl = [r for r in rows if r[1] != "self.link"]
    assert len(sl) == 77 and len(ctl) == 16
    sig = [r for r in sl if float(r[9]) < 0.05]
    assert len(sig) == 44
    assert sum(1 for r in sig if (float(r[2]) > 0) != (float(r[6]) > 0)) == 2   # "in only two of the 44 cases"
    assert sum(1 for r in ctl if float(r[9]) < 0.05) == 2


def test_ase_genes_are_self_linking_in_S1(S1, S2):
    sl = {r[0] for r in S1[1] if r[2] == "1"}
    assert {r[0] for r in S2[1] if r[1] == "self.link"} <= sl


def test_ase_genes_meet_selection(S1, S2):
    # chosen for p < 1e-8 and at least a 1.2-fold difference
    s1 = {r[0]: r for r in S1[1]}
    for r in S2[1]:
        if r[1] != "self.link":
            continue
        assert float(s1[r[0]][6]) > -math.log(1e-8) and abs(float(s1[r[0]][3])) >= math.log2(1.2), r[0]


@pytest.mark.xfail(strict=True, reason=(
    "Abstract: 'in 52%-78% of these genes, local regulatory variation acts directly in cis'. Table S2 "
    "(journal.pgen.0010025.st002.txt, column ase_pval) gives 44 of 77 self-linking genes with p < 0.05 (57%, as the "
    "Results say) and 42 of 77 (54.5%) when ASE also agrees in direction; no reading of the table gives 52%."))
def test_lower_bound_52_percent(S2):
    sl = [r for r in S2[1] if r[1] == "self.link"]
    assert round(100 * sum(1 for r in sl if float(r[9]) < 0.05 and (float(r[2]) > 0) == (float(r[6]) > 0)) / 77) == 52


def test_amn1_2_2_fold_by_higher(S1):
    r = [r for r in S1[1] if r[0] == "YBR158W"][0]
    assert r[2] == "1" and round(2 ** float(r[3]), 1) == 2.2


def test_lazy_table_matches_counts(L, S1):
    assert L["_name"] == "eqtl"
    n = len(L["id"])
    assert n == L["_n"]
    unplaced = 1428 - n
    assert 0 <= unplaced < 40
    assert all(v == 1 for v in L["loc"])


def test_positions_inside_R64_points(L, D):
    G = D["genes"]
    for i in range(len(L["id"])):
        c = L["c"][i]
        assert c in R64 and 1 <= L["p"][i] <= R64[c]
        assert L["a"][i] == L["e"][i] == L["p"][i]
        g = L["g"][i]
        assert G["c"][g] == c and L["p"][i] in (G["a"][g], G["e"][g])


def test_direction_follows_fold_change(L):
    for i in range(len(L["id"])):
        assert (L["up"][i] == "BY") == (L["ef"][i] > 0)


def test_dens_recount(D, ds, L):
    from collections import Counter
    hb = D["hb"]
    want = Counter((c, p // hb) for c, p in zip(L["c"], L["p"]))
    got = {(D["dens"]["c"][k], D["dens"]["b"][k]): D["dens"]["n"][k] for k in range(len(D["dens"]["ds"])) if D["dens"]["ds"][k] == ds}
    assert got == dict(want)
    assert {D["dens"]["k"][k] for k in range(len(D["dens"]["ds"])) if D["dens"]["ds"][k] == ds} == {"local"}


def test_provenance(D, ds, L):
    codes = {c for c, d in zip(D["sources"]["code"], D["sources"]["ds"]) if d == ds}
    assert codes == {"r05:S1", "r05:S2", "r05:TF"}
    for pv in L["pv"]:
        code, row = pv.split("#")
        assert code == "r05:S1" and int(row) >= 2


def test_verbatim_rows_kept(L):
    seen = 0
    for c in R64:
        with open(os.path.join(ROOT, "docs", "v", f"r05_{c}.json"), encoding="utf-8") as f:
            v = json.load(f)
        assert len(v["cols"]["r05:S1"]) == 9
        for rid, row in v["rows"]["r05:S1"].items():
            assert len(row) == 9
            seen += 1
    assert seen == len(L["id"])
    with open(os.path.join(ROOT, "docs", "v", "r05_extra.json"), encoding="utf-8") as f:
        x = json.load(f)
    assert len(x["r05:S1"]["rows"]) == 5727 and len(x["r05:S2"]["rows"]) == 93
