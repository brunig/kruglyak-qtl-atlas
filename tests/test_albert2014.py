"""Validate the ingested Albert et al. 2014 (Nature 506:494) X-pQTL against the paper and its files.

    python3 -I -m pytest -q tests/test_albert2014.py      (after ./build.sh)

The main text is paywalled; the paper's numbers come from the abstract, the Extended Data legends
and the Supplementary Information (quoted in sources/albert2014/extracted/text_facts.tsv).
Disagreements are strict xfails naming file and column; the data are not changed.
"""
import collections
import json
import os
import re

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "sources", "albert2014")
RAW = os.path.join(SRC, "raw")
P = "41586_2014_BFnature12904_"
R64 = {"chrI": 230218, "chrII": 813184, "chrIII": 316620, "chrIV": 1531933, "chrV": 576874, "chrVI": 270161,
       "chrVII": 1090940, "chrVIII": 562643, "chrIX": 439888, "chrX": 745751, "chrXI": 666816, "chrXII": 1078177,
       "chrXIII": 924431, "chrXIV": 784333, "chrXV": 1091291, "chrXVI": 948066}
ROMAN = list(R64)


@pytest.fixture(scope="module")
def D():
    with open(os.path.join(ROOT, "docs", "data.json"), encoding="utf-8") as f:
        return json.load(f)


@pytest.fixture(scope="module")
def ds(D):
    return D["datasets"]["id"].index("albert2014")


@pytest.fixture(scope="module")
def L(D, ds):
    k = [i for i, d in enumerate(D["lazy"]["ds"]) if d == ds and D["lazy"]["name"][i] == "pqtl"][0]
    with open(os.path.join(ROOT, "docs", D["lazy"]["file"][k]), encoding="utf-8") as f:
        t = json.load(f)
    t["_n"] = D["lazy"]["n"][k]
    return t


def _xlsx(name):
    import openpyxl
    wb = openpyxl.load_workbook(os.path.join(RAW, P + name), read_only=True, data_only=True)
    rows = [list(r) for r in wb.worksheets[0].iter_rows(values_only=True)]
    return rows[0], rows[1:]


@pytest.fixture(scope="module")
def SD1():
    h, rows = _xlsx("MOESM49_ESM.xlsx")
    return {c: i for i, c in enumerate(h)}, rows


@pytest.fixture(scope="module")
def SD2():
    h, rows = _xlsx("MOESM50_ESM.xlsx")
    return {c: i for i, c in enumerate(h)}, rows


@pytest.fixture(scope="module")
def ST():
    with open(os.path.join(SRC, "extracted", "supp_tables.tsv")) as f:
        lines = [l.rstrip("\n").split("\t") for l in f if not l.startswith("#")]
    head = lines[0]
    return [dict(zip(head, r)) for r in lines[1:]]


@pytest.fixture(scope="module")
def genes(D):
    return D["genes"]


def analysed(SD1):
    h, rows = SD1
    return {r[h["ORF"]] for r in rows if r[h["localExperiment?"]] in ("yes", "no")}


# ---------------------------------------------------------------- citation, design ---------
def test_citation_matches_article_page(D, ds):
    html = open(os.path.join(SRC, "article.html"), encoding="utf-8").read()
    authors = re.findall(r'<meta name="citation_author" content="([^"]*)"', html)
    assert authors == ["Albert, Frank W.", "Treusch, Sebastian", "Shockley, Arthur H.", "Bloom, Joshua S.", "Kruglyak, Leonid"]
    assert '<meta name="citation_volume" content="506"' in html and '<meta name="citation_firstpage" content="494"' in html
    assert D["datasets"]["doi"][ds] == "10.1038/nature12904"
    assert "Nature 506:494-497" in D["datasets"]["cite"][ds]


def test_quoted_facts_are_verbatim():
    import html as H
    raw = open(os.path.join(SRC, "article.html"), encoding="utf-8").read()
    raw = re.sub(r"<(script|style)[^>]*>.*?</\1>", " ", raw, flags=re.S)
    html = re.sub(r"\s+", " ", H.unescape(re.sub(r"<[^>]+>", " ", raw)))
    import pdfplumber
    si = re.sub(r"\s+", " ", "\n".join((p.extract_text() or "") for p in pdfplumber.open(os.path.join(RAW, P + "MOESM48_ESM.pdf")).pages))
    with open(os.path.join(SRC, "extracted", "text_facts.tsv")) as f:
        rows = [l.rstrip("\n").split("\t") for l in f][1:]
    assert len(rows) == 15
    for kind, where, quote, found in rows:
        q = re.sub(r"\s+", " ", quote)
        assert q in (si if found.endswith(".pdf") else html) or (q.replace("’", "'") in html), quote


def test_160_genes_analysed(SD1):
    # abstract: "We applied this method to 160 genes"
    h, rows = SD1
    assert len(rows) == 174
    st = collections.Counter(r[h["localExperiment?"]] for r in rows)
    assert st["yes"] + st["no"] == 160
    assert sum(v for k, v in st.items() if k not in ("yes", "no")) == 14


def test_41_local_experiments_match_supp_table_2(SD1, ST):
    h, rows = SD1
    local = {r[h["ORF"]] for r in rows if r[h["localExperiment?"]] == "yes"}
    st2 = {r["col1"] for r in ST if r["table"] == "ST2"}
    assert len(local) == 41 and st2 == local


def test_local_experiment_genes_named_in_extended_data_fig_4(SD1, genes):
    # "YMR315W and GCN1 were 'local' experiments"
    h, rows = SD1
    local = {r[h["ORF"]] for r in rows if r[h["localExperiment?"]] == "yes"}
    gcn1 = genes["id"][genes["n"].index("GCN1")]
    assert "YMR315W" in local and gcn1 in local


def test_leu1_has_14_x_pqtl(SD2, genes):
    # Extended Data Figure 3: "For this gene, 14 X-pQTL are called." (the gene is LEU1)
    h, rows = SD2
    leu1 = genes["id"][genes["n"].index("LEU1")]
    assert sum(1 for r in rows if r[h["gene"]] == leu1) == 14


def test_all_lod_above_threshold(L):
    # Extended Data Figure 3: genome-wide significance threshold lod = 4.5
    assert min(L["lod"]) >= 4.5


@pytest.mark.skip(reason="Supplementary Table 1 / Note 1 (27 replicated genes, 240 X-pQTL in detection sets) are computed from replicate count data in Supplementary Data 3; the atlas does not re-map")
def test_reproducibility_240():
    pass


def test_supp_note_2_701_is_consistent(SD1, SD2):
    """Supplementary Note 2: 701 X-pQTL that did not overlap a significant eQTL, were not on chrII or
    chrIII, and not on the gene's own chromosome. The eQTL overlap cannot be applied from the files
    (Supplementary Data 1 gives no eQTL positions), so only check that the other two filters leave at least 701."""
    h1, r1 = SD1
    own = {r[h1["ORF"]]: r[h1["chr"]] for r in r1}
    h, rows = SD2
    keep = analysed(SD1)
    n = sum(1 for r in rows if r[h["gene"]] in keep and r[h["chromosome"]] not in ("chr02", "chr03") and r[h["chromosome"]] != own[r[h["gene"]]])
    assert n >= 701


# ---------------------------------------------------------------- the lazy file -----------
def test_lazy_file_matches_table(L, SD2):
    _, rows = SD2
    assert L["_n"] == len(L["id"]) == len(rows) == 1025
    assert all(len(L[k]) == 1025 for k in ("g", "c", "a", "e", "p", "lod", "ef", "loc", "st", "pv"))


def test_coordinates_inside_r64(L):
    for c, a, e, p in zip(L["c"], L["a"], L["e"], L["p"]):
        assert c in R64 and 0 <= a <= p <= e <= R64[c]


def test_genes_resolve(L):
    assert min(L["g"]) >= 0


def test_effect_sign_as_given(L, SD2):
    h, rows = SD2
    for ef, r in zip(L["ef"], rows):
        assert ef == pytest.approx(r[h["alleleFrequencyDifference"]], rel=1e-5)


def test_status_flags_excluded_genes(L, SD1):
    keep = analysed(SD1)
    h, rows = SD1
    status = {r[h["ORF"]]: r[h["localExperiment?"]] for r in rows}
    gene = [None] * len(L["id"])
    import openpyxl  # noqa: F401
    h2, r2 = _xlsx("MOESM50_ESM.xlsx")
    for i, r in enumerate(r2):
        gene[i] = r[0]
    off = {g for g in gene if g not in keep}
    assert len(off) == 9 and sum(1 for g in gene if g not in keep) == 66
    assert all(L["st"][i] == ("analysed" if gene[i] in keep else status[gene[i]]) for i in range(len(gene)))


def test_local_by_overlap(L, SD1):
    h, rows = SD1
    pos = {r[h["ORF"]]: ("chr" + ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI"][int(r[h["chr"]][3:]) - 1],
                         min(r[h["startPosition"]], r[h["stopPosition"]]), max(r[h["startPosition"]], r[h["stopPosition"]])) for r in rows}
    _, r2 = _xlsx("MOESM50_ESM.xlsx")
    n = 0
    for i, r in enumerate(r2):
        c, lo, hi = pos[r[0]]
        want = 1 if (c == L["c"][i] and L["a"][i] <= hi and L["e"][i] >= lo) else 0
        assert L["loc"][i] == want
        n += want
    assert n == 5


def test_density_recounts_from_lazy_file(D, ds, L):
    hb = D["hb"]
    want = collections.Counter()
    for i in range(len(L["id"])):
        if L["st"][i] != "analysed":
            continue
        for b in range(L["a"][i] // hb, L["e"][i] // hb + 1):
            want[("local" if L["loc"][i] else "distant", L["c"][i], b)] += 1
    Dn = D["dens"]
    got = {(Dn["k"][k], Dn["c"][k], Dn["b"][k]): Dn["n"][k] for k in range(len(Dn["ds"])) if Dn["ds"][k] == ds}
    assert got == dict(want)


@pytest.mark.xfail(strict=True, reason=(
    "Supplementary Data 1 (MOESM49_ESM.xlsx) column number_of_trans_pQTL does not count this study's X-pQTL: it equals the "
    "number of rows per gene in Supplementary Data 2 (MOESM50_ESM.xlsx, column gene) for only 4 of the 91 analysed genes with a value. "
    "Most defensible reading: the column (like local_pQTL_LOD and the eQTL columns) reports earlier mass-spectrometry pQTL and eQTL studies used for comparison."))
def test_sd1_trans_pqtl_counts_equal_sd2(SD1, SD2):
    h1, r1 = SD1
    h2, r2 = SD2
    cnt = collections.Counter(r[h2["gene"]] for r in r2)
    keep = analysed(SD1)
    vals = [(r[h1["ORF"]], r[h1["number_of_trans_pQTL"]]) for r in r1 if r[h1["number_of_trans_pQTL"]] != "NA" and r[h1["ORF"]] in keep]
    assert all(cnt[g] == v for g, v in vals)


def test_sd1_trans_pqtl_counts_as_observed(SD1, SD2):
    h1, r1 = SD1
    h2, r2 = SD2
    cnt = collections.Counter(r[h2["gene"]] for r in r2)
    keep = analysed(SD1)
    vals = [(r[h1["ORF"]], r[h1["number_of_trans_pQTL"]]) for r in r1 if r[h1["number_of_trans_pQTL"]] != "NA" and r[h1["ORF"]] in keep]
    assert len(vals) == 91 and sum(1 for g, v in vals if cnt[g] == v) == 4


def test_supp_table_3_genes_have_a_chrII_x_pqtl(ST, SD2, genes):
    """Supplementary Table 3 lists 49 proteins affected by the chromosome II hotspot with an effect.
    Every one that matches an SGD standard name has a chrII X-pQTL in Supplementary Data 2 of the same sign;
    24 of the 48 equal the nearest chrII peak's allele frequency difference at 2 decimals (the table's effect is
    probably measured at the hotspot, not at each gene's own peak; recorded in OPEN_QUESTIONS)."""
    h, rows = SD2
    byname = {n: genes["id"][i] for i, n in enumerate(genes["n"]) if n}
    st3 = [(r["col1"], float(r["col2"])) for r in ST if r["table"] == "ST3"]
    assert len(st3) == 49
    unmatched, same2 = [], 0
    for name, eff in st3:
        orf = byname.get(name, name if re.match(r"^Y[A-P][LR]\d{3}[WC]", name) else None)
        cand = [r[h["alleleFrequencyDifference"]] for r in rows if r[h["gene"]] == orf and r[h["chromosome"]] == "chr02"]
        if not cand:
            unmatched.append(name)
            continue
        best = min(cand, key=lambda v: abs(round(v, 2) - eff))
        assert (best > 0) == (eff > 0), name
        same2 += round(best, 2) == eff
    assert unmatched == ["SUR4"] and same2 == 24


# ---------------------------------------------------------------- provenance -------------
def test_provenance_resolves(D, L):
    codes = set(D["sources"]["code"])
    assert all(pv.split("#")[0] in codes for pv in L["pv"])
    for pv in D["crosses"]["pv"]:
        assert pv.split("#")[0] in codes


def test_verbatim_rows_kept(L):
    seen = 0
    for c in R64:
        with open(os.path.join(ROOT, "docs", "v", f"a14_{c}.json"), encoding="utf-8") as f:
            v = json.load(f)
        assert len(v["cols"]["a14:SD2"]) == 7
        for pid, row in v["rows"]["a14:SD2"].items():
            assert len(row) == 7
            seen += 1
    assert seen == 1025
    with open(os.path.join(ROOT, "docs", "v", "a14_extra.json"), encoding="utf-8") as f:
        x = json.load(f)
    assert {k: len(v["rows"]) for k, v in x.items()} == {"a14:SD1": 174, "a14:ST": 95, "a14:TF": 15, "a14:SD3": 4}


def test_build_consistent_with_r64(SD1, genes):
    h, rows = SD1
    gi = {}
    for i, g in enumerate(genes["id"]):
        gi.setdefault(g, i)
    for r in rows:
        g = gi[r[h["ORF"]]]
        c = "chr" + ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI"][int(r[h["chr"]][3:]) - 1]
        lo, hi = sorted((r[h["startPosition"]], r[h["stopPosition"]]))
        assert genes["c"][g] == c and abs(genes["a"][g] - lo) <= 2000 and abs(genes["e"][g] - hi) <= 2000
