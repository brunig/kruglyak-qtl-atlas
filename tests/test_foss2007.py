"""Validate the ingested Foss et al. 2007 (Nat Genet 39:1369) tables against the paper's supplement.

    python3 -I -m pytest -q tests/test_foss2007.py      (after python3 -I ingest/foss2007.py)
"""
import json
import os

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "sources", "foss2007")
R64 = {"chrI": 230218, "chrII": 813184, "chrIII": 316620, "chrIV": 1531933, "chrV": 576874, "chrVI": 270161,
       "chrVII": 1090940, "chrVIII": 562643, "chrIX": 439888, "chrX": 745751, "chrXI": 666816, "chrXII": 1078177,
       "chrXIII": 924431, "chrXIV": 784333, "chrXV": 1091291, "chrXVI": 948066}


@pytest.fixture(scope="module")
def D():
    with open(os.path.join(ROOT, "docs", "data.json"), encoding="utf-8") as f:
        return json.load(f)


@pytest.fixture(scope="module")
def ds(D):
    return D["datasets"]["id"].index("foss2007")


@pytest.fixture(scope="module")
def L():
    with open(os.path.join(ROOT, "docs", "v", "f07_mspqtl.json"), encoding="utf-8") as f:
        return json.load(f)


def tsv(name):
    with open(os.path.join(SRC, "extracted", name), encoding="utf-8") as f:
        lines = f.read().splitlines()
    head = lines[1].split("\t")
    return [dict(zip(head, l.split("\t"))) for l in lines[2:]]


@pytest.fixture(scope="module")
def facts():
    return {r["key"]: r["quote"] for r in tsv("text_facts.tsv")}


def ints(s):
    import re
    return [int(x) for x in re.findall(r"\d+", s)]


def test_citation(D, ds):
    assert D["datasets"]["doi"][ds] == "10.1038/ng.2007.22"
    assert "Nat Genet 39:1369" in D["datasets"]["cite"][ds]


def test_quotes_verbatim_in_sources():
    import html
    import re
    import pdfplumber
    with pdfplumber.open(os.path.join(SRC, "raw", "41588_2007_BFng200722_MOESM17_ESM.pdf")) as pdf:
        si = re.sub(r"\s+", " ", " ".join((p.extract_text() or "") for p in pdf.pages))
    art = re.sub(r"\s+", " ", html.unescape(re.sub(r"<[^>]+>", " ", open(os.path.join(SRC, "article.html"), encoding="utf-8", errors="replace").read())))
    for r in tsv("text_facts.tsv"):
        assert r["quote"] in (si if r["source"] == "si" else art), r["key"]


def test_one_cross_by_rm(D, ds):
    X = [i for i, d in enumerate(D["crosses"]["ds"]) if d == ds]
    assert len(X) == 1
    P = D["parents"]
    assert {P["name"][D["crosses"]["p1"][X[0]]], P["name"][D["crosses"]["p2"][X[0]]]} == {"BY", "RM"}


def test_six_local_linkages(L, facts):
    # "6 peptides showed linkage to a marker that was less than 20 kilobases from the gene encoding that peptide."
    assert "6 peptides showed linkage" in facts["st1_caption"]
    assert len(L["id"]) == 6 and all(v == 1 for v in L["loc"])
    assert all(d < 20000 for d in L["dist"])


def test_lods_above_three(L):
    assert all(l > 3 for l in L["lod"])


def test_local_marker_on_gene_chromosome(L, D):
    """ORF-named local markers lie on the gene's chromosome in R64-5-1 and inside the placed interval;
    chromosome-arm markers (gXLnn) name the gene's chromosome by letter."""
    G = D["genes"]
    roman = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI"]
    checked = 0
    for i, m in enumerate(L["mk"]):
        name = m.rsplit("_", 1)[0]
        if name in G["id"]:
            k = G["id"].index(name)
            assert G["c"][k] == L["c"][i]
            assert L["a"][i] <= (G["a"][k] + G["e"][k]) // 2 <= L["e"][i] + 2000
            checked += 1
        elif name.startswith("g") and name[2] in "LR":
            assert "chr" + roman[ord(name[1]) - 65] == L["c"][i]
            checked += 1
    assert checked == 5   # the sixth (NAL002C_307) is neither an SGD name nor an arm marker: OPEN_QUESTIONS


def test_placement_inside_r64(L):
    for i in range(len(L["id"])):
        assert 1 <= L["a"][i] <= L["p"][i] <= L["e"][i] <= R64[L["c"][i]]


def test_supplementary_note_sums(facts):
    p = ints(facts["peptide_links"]); assert (p[0], p[1], p[3]) == (71, 163, 93)       # "... LOD > 3, ... total of 93"
    q = ints(facts["transcript_links"]); assert (q[0], q[1], q[3]) == (110, 163, 140)
    o = ints(facts["peptide_overlap"])
    assert o[0] == 93 and o[1] + o[3] + o[6] == 93          # 18 + 19 + 56
    t = ints(facts["transcript_overlap"])
    assert t[0] == 140 and t[1] + t[3] + t[6] == 140        # 17 + 31 + 92
    assert round(100 * 59 / 137) == 43


def test_no_hotspot_rows(D, ds):
    assert sum(1 for d in D["qtl"]["ds"] if d == ds) == 0


def test_lazy_row_and_density(D, ds, L):
    k = [i for i, d in enumerate(D["lazy"]["ds"]) if d == ds]
    assert len(k) == 1 and D["lazy"]["name"][k[0]] == "mspqtl" and D["lazy"]["n"][k[0]] == len(L["id"])
    from collections import Counter
    rec = Counter()
    for i in range(len(L["id"])):
        for b in range(L["a"][i] // D["hb"], L["e"][i] // D["hb"] + 1):
            rec[(L["c"][i], b)] += 1
    got = {(D["dens"]["c"][j], D["dens"]["b"][j]): D["dens"]["n"][j] for j in range(len(D["dens"]["ds"])) if D["dens"]["ds"][j] == ds}
    assert got == dict(rec)


def test_provenance_and_verbatim_rows(D, ds, L):
    codes = {c for c, d in zip(D["sources"]["code"], D["sources"]["ds"]) if d == ds}
    for pv in L["pv"]:
        c, r = pv.split("#")
        assert c in codes and int(r) >= 2
    seen = 0
    for ch in R64:
        with open(os.path.join(ROOT, "docs", "v", f"f07_{ch}.json"), encoding="utf-8") as f:
            v = json.load(f)
        seen += len(v["rows"]["f07:ST1"])
        assert all(len(r) == len(v["cols"]["f07:ST1"]) for r in v["rows"]["f07:ST1"].values())
    assert seen == 6
