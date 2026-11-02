"""Validate the ingested Gatbonton et al. 2006 tables against the paper and its supplement.

    python3 -I -m pytest -q tests/test_gatbonton2006.py      (after the ingest)

Claims from the paper are asserted as stated; where the files disagree, the test is a strict xfail
naming the file and column.
"""
import collections
import html
import json
import os
import re

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "sources", "gatbonton2006")
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
    return D["datasets"]["id"].index("gatbonton2006")


@pytest.fixture(scope="module")
def Q(D, ds):
    return [r for r in rows(D["qtl"]) if r["ds"] == ds]


@pytest.fixture(scope="module")
def CG(D, ds):
    return [r for r in rows(D["causal"]) if r["ds"] == ds]


@pytest.fixture(scope="module")
def facts():
    with open(os.path.join(SRC, "extracted", "text_facts.tsv")) as f:
        r = [l.rstrip("\n").split("\t") for l in f]
    return [dict(zip(r[0], x)) for x in r[1:]]


@pytest.fixture(scope="module")
def S1():
    with open(os.path.join(SRC, "extracted", "table_s1.tsv")) as f:
        r = [l.rstrip("\n").split("\t") for l in f if not l.startswith("#")]
    return r[1:]


def page(name):
    t = open(os.path.join(SRC, name), encoding="utf-8", errors="replace").read()
    m = re.search(r'<div class="article-text"[^>]*>(.*)', t, re.S)
    t = m.group(1) if m else t
    t = re.sub(r"<(script|style)[^>]*>.*?</\1>", " ", t, flags=re.S)
    return re.sub(r"\s+", " ", html.unescape(re.sub(r"<[^>]+>", " ", t)))


def genes(D):
    G = D["genes"]
    return {G["n"][i]: (G["c"][i], G["a"][i]) for i in range(len(G["id"])) if G["n"][i]} | \
           {G["id"][i]: (G["c"][i], G["a"][i]) for i in range(len(G["id"]))}


def test_quotes_verbatim(facts):
    pages = {"article.html": page("article.html"), "correction.html": page("correction.html")}
    assert len(facts) == 20
    for f in facts:
        assert f["quote"] in pages[f["source"]], f["quote"][:60]


def test_citation(D, ds):
    assert D["datasets"]["doi"][ds] == "10.1371/journal.pgen.0020035"
    t = page("article.html")
    assert "Telomere Length as a Quantitative Trait" in open(os.path.join(SRC, "article.html"), encoding="utf-8").read()
    assert "PLoS Genet 2(6): e104" in page("correction.html")


def test_one_cross_122_segregants(D, ds):
    X = [r for r in rows(D["crosses"]) if r["ds"] == ds]
    assert len(X) == 1 and X[0]["n"] == 122
    P = D["parents"]
    assert {P["name"][X[0]["p1"]], P["name"][X[0]["p2"]]} == {"BY", "RM"}


def test_two_qtl_chr12_chr13(Q):
    assert [(r["c"], r["lod"], r["ve"]) for r in Q] == [("chrXII", 7.4, 0.25), ("chrXIII", 4.3, 0.13)]


def test_lod_above_permutation_max(Q):
    # "we did not observe any scores higher than four in the analysis of 100 permuted datasets"
    assert all(r["lod"] > 4 for r in Q)


def test_chromosome_only(Q):
    assert all(r["a"] is None and r["e"] is None and r["p"] is None and r["c"] in R64 for r in Q)


def test_direction(D, Q):
    P = D["parents"]["name"]
    assert P[Q[0]["up"]] == "RM"      # chr12: RM11 allele confers longer telomeres
    assert Q[1]["up"] == -1           # chr13: direction only in the Table 2 image


def test_candidates_named_in_text_on_qtl_chromosomes(D, CG, Q):
    assert sorted(r["name"] for r in CG) == sorted(["EST1", "VPS34", "RPS28B", "ARV1", "VPS9"])
    ids = {i: r for i, r in enumerate(rows(D["qtl"]))}
    for r in CG:
        assert r["tier"] == "c" and r["c"] == ids[r["q1"]]["c"]


def test_s1_genes_with_substitutions(S1):
    aa = {r[1].strip(): r[4].strip() for r in S1 if len(r) > 4}
    assert "11 AA substitutions" in aa["EST1"] and "7 AA substitutions" in aa["VPS9"]


@pytest.mark.xfail(strict=True, reason=(
    "Table S1 (journal.pgen.0020035.st001.xls, Sheet1, column B) lists BET1 under 'CHROM 13', but BET1 is YIL004C on "
    "chrIX in SGD R64-5-1. The other 13 genes lie on the chromosome they are listed under."))
def test_s1_genes_on_listed_chromosome(D, S1):
    G = genes(D)
    chrom = None
    bad = []
    for r in S1:
        name = r[1].strip().replace(" ", "")
        if name.startswith("CHROM"):
            chrom = {"12": "chrXII", "13": "chrXIII"}[name[5:]]
            continue
        if chrom and name in G and G[name][0] != chrom:
            bad.append(name)
    assert not bad


def test_s1_other_genes_on_listed_chromosome(D, S1):
    G = genes(D)
    chrom, n = None, 0
    for r in S1:
        name = r[1].strip().replace(" ", "")
        if name.startswith("CHROM"):
            chrom = {"12": "chrXII", "13": "chrXIII"}[name[5:]]
            continue
        if chrom and name in G and name != "BET1":
            assert G[name][0] == chrom, name
            n += 1
    assert n == 13


def test_survey_counts_in_text(facts):
    q = " ".join(f["quote"] for f in facts if f["kind"] == "survey")
    assert "72 genes" in q and "80 genes" in q and "88 have not" in q and "39 genes" in q


def test_correction_recorded(facts):
    assert any(f["kind"] == "correction" and "Table 1 was incomplete" in f["quote"] for f in facts)


def test_provenance(D, ds):
    codes = {c for c, d in zip(D["sources"]["code"], D["sources"]["ds"]) if d == ds}
    for t in ("crosses", "qtl", "causal"):
        for pv, d in zip(D[t]["pv"], D[t]["ds"]):
            if d != ds:
                continue
            code, row = pv.split("#")
            assert code in codes and int(row) >= 2


def test_original_rows_kept(D, ds):
    ids = {r["id"] for r in rows(D["qtl"]) if r["ds"] == ds} | {r["id"] for r in rows(D["causal"]) if r["ds"] == ds}
    seen = set()
    for ch in R64:
        with open(os.path.join(ROOT, "docs", "v", f"g06_{ch}.json"), encoding="utf-8") as f:
            v = json.load(f)
        seen |= set(v["rows"].get("g06:TF", {}))
    assert seen == ids
    with open(os.path.join(ROOT, "docs", "v", "g06_extra.json"), encoding="utf-8") as f:
        ex = json.load(f)
    assert len(ex["g06:TF"]["rows"]) == 20 and len(ex["g06:S1"]["rows"]) == 17
