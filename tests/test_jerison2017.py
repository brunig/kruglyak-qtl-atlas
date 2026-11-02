"""Validate the ingested Jerison et al. 2017 tables against the paper and its supplement.

    python3 -I -m pytest -q tests/test_jerison2017.py      (after the ingest)

Claims come from the PMC full text (sources/jerison2017/article_pmc5580887.xml) as quoted in
sources/jerison2017/extracted/text_facts.tsv. Where the supplement files disagree with the text,
the test is a strict xfail naming the file and column.
"""
import csv
import json
import os
import re
import xml.etree.ElementTree as ET

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "sources", "jerison2017")
RAW = os.path.join(SRC, "raw")
P = "elife_poa_e27167_"
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
    return D["datasets"]["id"].index("jerison2017")


@pytest.fixture(scope="module")
def Q(D, ds):
    return [r for r in rows(D["qtl"]) if r["ds"] == ds]


@pytest.fixture(scope="module")
def text():
    r = ET.parse(os.path.join(SRC, "article_pmc5580887.xml")).getroot()
    norm = lambda n: " ".join(" ".join(n.itertext()).split()) if n is not None else ""
    return norm(r.find(".//abstract")) + " " + norm(r.find(".//body")) + " " + norm(r.find(".//back"))


@pytest.fixture(scope="module")
def facts():
    with open(os.path.join(SRC, "extracted", "text_facts.tsv")) as f:
        lines = [l.rstrip("\n").split("\t") for l in f if not l.startswith("#")]
    return {k: q for k, _, q in lines[1:]}


@pytest.fixture(scope="module")
def sf5():
    with open(os.path.join(RAW, P + "Supplementary_File_5.csv"), newline="") as f:
        return list(csv.reader(f))


@pytest.fixture(scope="module")
def sf2():
    with open(os.path.join(RAW, P + "Supplementary_File_2.csv"), newline="") as f:
        return list(csv.reader(f, delimiter=";"))


@pytest.fixture(scope="module")
def sf3():
    with open(os.path.join(RAW, P + "Supplementary_File_3.csv"), newline="") as f:
        return list(csv.DictReader(f))


def test_citation_matches_article(D, ds):
    r = ET.parse(os.path.join(SRC, "article_pmc5580887.xml")).getroot()
    title = "".join(r.find(".//article-title").itertext())
    assert title == "Genetic variation in adaptability and pleiotropy in budding yeast"
    assert r.findtext(".//volume") == "6" and r.findtext(".//elocation-id") == "e27167"
    authors = [n.findtext("surname") for n in r.findall('.//contrib[@contrib-type="author"]/name')]
    assert authors == ["Jerison", "Kryazhimskiy", "Mitchell", "Bloom", "Kruglyak", "Desai"]
    assert D["datasets"]["doi"][ds] == "10.7554/eLife.27167"


def test_quotes_are_verbatim(text, facts):
    for k, q in facts.items():
        assert q in text, k


def test_230_founders(D, ds, sf2, facts):
    assert "230 segregants" in facts["founders"]
    assert len(sf2) - 1 == 230
    with open(os.path.join(RAW, P + "Supplementary_File_1.csv")) as f:
        assert sum(1 for _ in f) - 1 == 230
    x = [r for r in rows(D["crosses"]) if r["ds"] == ds]
    assert len(x) == 1 and x[0]["n"] == 230


def test_1840_populations(facts):
    assert "1840" in facts["populations"] and 230 * 8 == 1840


def _lines(sf2):
    ot = sum(len([x for x in r[5].split(",") if x.strip()]) for r in sf2[1:])
    ht = sum(len([x for x in r[6].split(",") if x.strip()]) for r in sf2[1:])
    return ot, ht


@pytest.mark.xfail(strict=True, reason=(
    "Text (Methods, Experimental evolution): 910 evolved lines at OT and 839 at HT. Supplementary file 2 "
    "(columns '(final fitness at 30C, final fitness at 37C), populations evolved at 30C/37C') lists 907 and 852 "
    "evolved populations."))
def test_910_and_839_lines(sf2):
    assert _lines(sf2) == (910, 839)


def test_lines_as_observed(sf2):
    assert _lines(sf2) == (907, 852)


def test_8_loci(sf5, facts, Q):
    assert "8 QTLs" in facts["n_qtl"]
    assert len(sf5) - 1 == 8
    assert len({q["lv"] for q in Q}) == 8


def test_qtl_rows_equal_y_cells(sf5, Q):
    head = sf5[0]
    ys = sum(1 for r in sf5[1:] for j, h in enumerate(head) if h.startswith("QTL for trait") and r[j].strip().lower() == "y")
    assert ys == len(Q) == 30


@pytest.mark.xfail(strict=True, reason=(
    "Text (Methods): 'the list of 99 genes within QTL confidence intervals' (dubious ORFs excluded). Supplementary "
    "file 5, column 'Genes within confidence intervals', lists 208 unique gene names over the 8 loci."))
def test_99_genes_in_intervals(sf5):
    genes = {g.strip() for r in sf5[1:] for g in re.split(r"[,\s]+", r[1]) if g.strip()}
    assert len(genes) == 99


def test_genes_in_intervals_as_observed(sf5):
    genes = {g.strip() for r in sf5[1:] for g in re.split(r"[,\s]+", r[1]) if g.strip()}
    assert len(genes) == 208


@pytest.mark.xfail(strict=True, reason=(
    "Text: 603 mutations fixed or above 50% in 254 sequenced lines. Supplementary file 3 ('All mutations called at "
    "frequencies of >= 50%') has 738 rows from 219 populations (column 'Population')."))
def test_603_mutations_254_lines(sf3):
    assert len(sf3) == 603 and len({r["Population"] for r in sf3}) == 254


def test_mutations_as_observed(sf3):
    assert len(sf3) == 738 and len({r["Population"] for r in sf3}) == 219
    assert all(float(r["Estimated frequency within population"]) >= 0.5 for r in sf3)


def test_coordinates_inside_R64(Q):
    for q in Q:
        assert q["c"] in R64 and 1 <= q["a"] <= q["p"] <= q["e"] <= R64[q["c"]]
        assert q["a"] == q["p"] == q["e"]    # no interval bounds are tabulated


def test_genes_listed_under_each_locus_lie_near_it(D, sf5):
    """Build check: the genes Supplementary file 5 lists under a locus sit on its chromosome, near it, on R64-5-1."""
    G = D["genes"]
    pos = {}
    for i, n in enumerate(G["n"]):
        if n:
            pos.setdefault(n, (G["c"][i], (G["a"][i] + G["e"][i]) // 2))
    for i, s in enumerate(G["id"]):
        pos.setdefault(s, (G["c"][i], (G["a"][i] + G["e"][i]) // 2))
    far, unresolved, n = [], [], 0
    for r in sf5[1:]:
        m = re.match(r"Chr([IVX]+) (\d+)", r[0]); c, p = "chr" + m.group(1), int(m.group(2))
        for g in {g.strip() for g in re.split(r"[,\s]+", r[1]) if g.strip()}:
            if g not in pos:
                unresolved.append(g); continue
            n += 1
            if not (pos[g][0] == c and abs(pos[g][1] - p) <= 60000):
                far.append((r[0], g))
    # Recorded naming issues in the source (OPEN_QUESTIONS Q4): MHF2 and GSY1 resolve to other chromosomes
    # under R64-5-1 names, MCP2 is not a current SGD name. Every other listed gene is within 60 kb.
    assert sorted(far) == [("ChrXII 646707", "GSY1"), ("ChrXV 154799", "MHF2")]
    assert unresolved == ["MCP2"] and n == 207


def test_kre33_candidate(D, ds, Q):
    G = D["genes"]; g = G["n"].index("KRE33")
    kq = [q for q in Q if q["lv"] == "ChrXIV 376315"]
    assert kq and all(q["g"] == g and q["gt"] == "c" for q in kq)
    assert G["c"][g] == "chrXIV" and abs((G["a"][g] + G["e"][g]) / 2 - 376315) < 10000
    assert sum(1 for q in Q if q["g"] >= 0) == len(kq)


def test_sign_reading_kre33_positive_in_adaptability_models(D, ds, Q):
    T = D["traits"]
    for name in ("Adaptability (fitness gain), OT", "Adaptability (fitness gain), HT"):
        t = T["name"].index(name)
        k = [q for q in Q if q["t"] == t and q["lv"] == "ChrXIV 376315"][0]
        assert k["ef"] > 0 and D["parents"]["name"][k["up"]] == "RM"


def test_coefficients(D, Q):
    """19 coefficients attached (the ingest refuses any model-table term on a cell not flagged 'y' in SF5);
    initial-fitness QTL have none in the bundle."""
    assert sum(1 for q in Q if q["ef"] is not None) == 19
    assert all(q["ef"] is None for q in Q if D["traits"]["name"][q["t"]].startswith("Initial fitness"))


def test_locus_matches_recorded():
    with open(os.path.join(ROOT, "docs", "v", "j17_extra.json")) as f:
        m = json.load(f)["j17:locus_matches"]["rows"]
    assert sorted(m) == [["ChrIV 1091232", "ChrIV 1062759", 28473], ["ChrXII 644077", "ChrXII 646707", 2630], ["ChrXIV 466590", "ChrXIV 470303", 3713]]


def test_table1_heritability_matches_text(D, ds):
    V = [r for r in rows(D["vexp"]) if D["traits"]["ds"][r["t"]] == ds]
    names = {D["traits"]["name"][r["t"]]: r for r in V}
    assert round(names["Adaptability (fitness gain), OT"]["h2"], 2) == 0.46
    assert round(names["Adaptability (fitness gain), HT"]["h2"], 2) == 0.44
    assert len(V) == 6


def test_provenance_resolves(D, ds):
    codes = {c for c, d in zip(D["sources"]["code"], D["sources"]["ds"]) if d == ds}
    for t in ("qtl", "causal", "vexp", "crosses"):
        for i, pv in enumerate(D[t]["pv"]):
            dsv = D[t]["ds"][i] if "ds" in D[t] else D["traits"]["ds"][D[t]["t"][i]]
            if dsv != ds:
                continue
            code, row = pv.split("#")
            assert code in codes and int(row) >= 1


def test_verbatim_rows_kept(Q):
    seen = set()
    for ch in R64:
        with open(os.path.join(ROOT, "docs", "v", f"j17_{ch}.json")) as f:
            v = json.load(f)
        for qid, row in v["rows"].get("j17:SF5", {}).items():
            assert len(row) == 10
            seen.add(qid)
    assert seen == {q["id"] for q in Q}
