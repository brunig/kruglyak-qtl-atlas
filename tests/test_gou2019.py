"""Validate the ingested Gou et al. 2019 tables against the paper (Genetics 211:731-740).

    python3 -I -m pytest -q tests/test_gou2019.py      (after ./build.sh or python3 -I ingest/gou2019.py)

Claims from the paper are asserted as stated; where the files disagree, the test is a strict
xfail naming file and column. Quotes in sources/gou2019/extracted/text_facts.tsv are checked
verbatim against the article PDF (per column, whitespace and line-break hyphens normalized).
"""
import collections
import json
import os
import re
import statistics

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "sources", "gou2019")
RAW = os.path.join(SRC, "raw")
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
    return D["datasets"]["id"].index("gou2019")


@pytest.fixture(scope="module")
def Q(D, ds):
    return [r for r in rows(D["qtl"]) if r["ds"] == ds]


@pytest.fixture(scope="module")
def T(D):
    return D["traits"]


def norm(s):
    return re.sub(r"\s+", "", s.replace("-\n", ""))


@pytest.fixture(scope="module")
def article():
    pdfplumber = pytest.importorskip("pdfplumber")
    cols, full = [], []
    with pdfplumber.open(os.path.join(RAW, "Gou_et_al_2019_Genetics.pdf")) as pdf:
        for pg in pdf.pages:
            full.append(pg.extract_text(x_tolerance=1) or "")
            for x0, x1 in ((0, pg.width / 2), (pg.width / 2, pg.width)):
                cols.append(pg.crop((x0, 0, x1, pg.height)).extract_text(x_tolerance=1) or "")
    return {"cols": norm("\n".join(cols)), "full": norm("\n".join(full)), "page1": full[0]}


@pytest.fixture(scope="module")
def supp():
    pdfplumber = pytest.importorskip("pdfplumber")
    with pdfplumber.open(os.path.join(RAW, "Supplemental Figures and Tables.pdf")) as pdf:
        return "\n".join(pg.extract_text() or "" for pg in pdf.pages)


def facts():
    with open(os.path.join(SRC, "extracted", "text_facts.tsv")) as f:
        r = [l.rstrip("\n").split("\t") for l in f]
    return [dict(zip(r[0], x)) for x in r[1:]]


def s2_means():
    import openpyxl
    wb = openpyxl.load_workbook(os.path.join(RAW, "S2.File_mutation_rate_engineer_strains.xlsx"), read_only=True, data_only=True)
    d = collections.defaultdict(list)
    for r in list(wb.worksheets[0].iter_rows(values_only=True))[1:]:
        d[r[1].strip()].append(r[2])
    return d


# ------------------------------------------------------------------ citation and design
def test_citation_from_pdf(D, ds, article):
    p1 = article["page1"]
    assert "doi: https://doi.org/10.1534/genetics.118.301609" in p1
    assert "Genetics, Vol. 211, 731–740 February 2019" in p1
    assert D["datasets"]["doi"][ds] == "10.1534/genetics.118.301609"
    assert "211:731-740" in D["datasets"]["cite"][ds]


def test_quotes_are_verbatim(article):
    nb = json.load(open(os.path.join(SRC, "code", "mutation rate qtl mapping.ipynb")))
    code = norm("".join("".join(c["source"]) for c in nb["cells"]))
    for f in facts():
        q = norm(f["quote"])
        assert q in article["cols"] or q in article["full"] or q in code, f["quote"]


def test_one_cross_843_of_1040_segregants(D, ds):
    X = [r for r in rows(D["crosses"]) if r["ds"] == ds]
    assert len(X) == 1 and X[0]["n"] == 843
    P = D["parents"]
    assert {P["name"][X[0]["p1"]], P["name"][X[0]["p2"]]} == {"BY", "RM"}
    assert 1040 - 843 == 197   # "we removed these 197 segregants"


def test_file_s1_has_843_scores():
    import openpyxl
    wb = openpyxl.load_workbook(os.path.join(RAW, "S1.File_mutation_rate_score.xlsx"), read_only=True, data_only=True)
    r = list(wb.worksheets[0].iter_rows(values_only=True))
    assert len(r) - 1 == 843 and all(isinstance(v[1], float) for v in r[1:])


def test_table_s2_groups_sum_to_1040(supp):
    m = re.search(r"Number of segregants (\d+) (\d+) (\d+) (\d+)", supp)
    assert sum(int(v) for v in m.groups()) == 1040


# ------------------------------------------------------------------ QTL
def mrate(Q, T):
    return {r["c"]: r for r in Q if T["raw"][r["t"]] == "mutation rate score"}


def test_four_mutation_rate_qtl(Q, T):
    assert set(mrate(Q, T)) == {"chrXII", "chrXIV", "chrI", "chrV"}


def test_variance_explained(Q, T):
    m = mrate(Q, T)
    assert [m[c]["ve"] for c in ("chrXII", "chrXIV", "chrI", "chrV")] == [0.088, 0.061, 0.031, 0.026]
    assert sum(r["ve"] for r in m.values()) == pytest.approx(0.207, abs=0.0015)   # text: 20.7% (parts sum to 20.6, rounding)
    assert 0.207 / 0.30 == pytest.approx(0.69, abs=0.005)                           # "69% of the estimated additive heritability"


def test_direction(Q, T, D):
    m, P = mrate(Q, T), D["parents"]
    assert {c: P["name"][r["up"]] for c, r in m.items()} == {"chrXII": "RM", "chrI": "RM", "chrXIV": "BY", "chrV": "BY"}


def test_rad5_mkt1_confirmed(Q, T, D):
    m, G = mrate(Q, T), D["genes"]
    for c, gene in (("chrXII", "RAD5"), ("chrXIV", "MKT1")):
        g = m[c]["g"]
        assert G["n"][g] == gene and G["c"][g] == c and m[c]["gt"] == "x"
    assert m["chrI"]["g"] == -1 and m["chrV"]["g"] == -1


def test_mutagen_qtl_overlap_xii_xiv(Q, T):
    for name in ("4NQO;0.1 mg/ml", "MMS;0.01%"):
        assert sorted(r["c"] for r in Q if T["raw"][r["t"]] == name) == ["chrXII", "chrXIV"]
    assert not [r for r in Q if T["raw"][r["t"]] == "H2O2;4 mM"]   # H2O2 loci are not identified in the text


def test_no_coordinates_anywhere(Q):
    """The files give chromosomes only; nothing is placed."""
    assert all(r["a"] is None and r["e"] is None and r["p"] is None and r["c"] in R64 for r in Q)


def test_thresholds_in_method(D, ds):
    for v in ("3.52", "3.62", "3.61", "3.64"):
        assert v in D["datasets"]["method"][ds]


# ------------------------------------------------------------------ engineered strains (File S2)
def test_file_s2_six_strains_96_each():
    d = s2_means()
    assert sorted(d) == sorted(["BY", "RM", "BY:RAD5", "RM:MKT1", "BY:RAD5-I791S", "BY:RAD5-E783D"])
    assert all(len(v) == 96 for v in d.values())


@pytest.mark.parametrize("a,b,stated", [("BY:RAD5-I791S", "BY", 0.9), ("BY:RAD5-E783D", "BY", 0.3), ("BY:RAD5-I791S", "BY:RAD5-E783D", 0.6)])
def test_variant_mean_differences(a, b, stated):
    d = s2_means()
    assert round((statistics.mean(d[a]) - statistics.mean(d[b])) / 1e-7, 1) == stated


@pytest.mark.xfail(strict=True, reason=(
    "Paper: BY::RAD5-RM minus BY = 2.9 x 10^-7. File S2 (S2.File_mutation_rate_engineer_strains.xlsx, columns "
    "'strain ' and 'mutation rate'): 4.62e-7 - 1.78e-7 = 2.84e-7, which rounds to 2.8."))
def test_rad5_replacement_difference():
    d = s2_means()
    assert round((statistics.mean(d["BY:RAD5"]) - statistics.mean(d["BY"])) / 1e-7, 1) == 2.9


@pytest.mark.xfail(strict=True, reason=(
    "Paper: RM::MKT1-BY minus RM = 6.1 x 10^-7. File S2: 8.27e-7 - 2.10e-7 = 6.17e-7, which rounds to 6.2."))
def test_mkt1_replacement_difference():
    d = s2_means()
    assert round((statistics.mean(d["RM:MKT1"]) - statistics.mean(d["RM"])) / 1e-7, 1) == 6.1


def test_rm_rate_differs_between_table_s1_and_file_s2(supp):
    """Recorded, not a test of the paper's claim: Table S1 gives RM 5.8e-7; File S2's RM replicates average 2.1e-7."""
    assert "RM11-1a 5.8×10-7" in supp
    assert round(statistics.mean(s2_means()["RM"]) / 1e-7, 1) == 2.1


# ------------------------------------------------------------------ provenance
def test_provenance_resolves(D, ds):
    codes = {c for c, d in zip(D["sources"]["code"], D["sources"]["ds"]) if d == ds}
    for t in ("qtl", "crosses", "causal", "vexp"):
        for i, pv in enumerate(D[t]["pv"]):
            if pv and pv.startswith("g19:"):
                code, row = pv.split("#")
                assert code in codes and int(row) >= 2
    assert len({c for c in D["sources"]["code"]}) == len(D["sources"]["code"])


def test_source_rows_kept(D, ds, Q):
    ids = {r["id"] for r in Q} | {r["id"] for r in rows(D["causal"]) if r["ds"] == ds}
    seen = set()
    for ch in R64:
        with open(os.path.join(ROOT, "docs", "v", f"g19_{ch}.json"), encoding="utf-8") as f:
            v = json.load(f)
        for rid, val in v["rows"]["g19:TF"].items():
            assert len(val) == len(v["cols"]["g19:TF"])
            seen.add(rid)
    assert seen == ids
    with open(os.path.join(ROOT, "docs", "v", "g19_extra.json"), encoding="utf-8") as f:
        x = json.load(f)
    assert len(x["g19:S1"]["rows"]) == 843 and len(x["g19:S2"]["rows"]) == 576 and len(x["g19:TF"]["rows"]) == len(facts())
