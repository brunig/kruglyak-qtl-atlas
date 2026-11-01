"""Validate the ingested Breunig et al. 2014 tables against the paper and against themselves.

    python3 -I -m pytest -q tests/test_breunig2014.py      (after ./build.sh)

Numbers come from the article text (sources/breunig2014/article.html) and Table S1's caption.
Disagreements are strict xfails naming the file and column; the data are not changed.
"""
import collections
import json
import os

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
S1 = "sources/breunig2014/raw/journal.pgen.1004142.s006.pdf, Table S1"


def rows(tab):
    keys = list(tab)
    return [dict(zip(keys, v)) for v in zip(*(tab[k] for k in keys))]


@pytest.fixture(scope="module")
def D():
    with open(os.path.join(ROOT, "docs", "data.json"), encoding="utf-8") as f:
        return json.load(f)


@pytest.fixture(scope="module")
def ds(D):
    return D["datasets"]["id"].index("breunig2014")


@pytest.fixture(scope="module")
def Q(D, ds):
    return [r for r in rows(D["qtl"]) if r["ds"] == ds]


@pytest.fixture(scope="module")
def T(D, ds):
    return {i: r for i, r in enumerate(rows(D["traits"])) if r["ds"] == ds}


@pytest.fixture(scope="module")
def raw(D, ds):
    """Table S1 rows verbatim, from the lazy files."""
    out = {}
    for ch in D["chroms"]:
        with open(os.path.join(ROOT, "docs", "v", f"b14_{ch}.json"), encoding="utf-8") as f:
            out.update(json.load(f)["rows"]["b14:S1"])
    return out


def per_compound(Q):
    return collections.Counter(r["t"] for r in Q)


# ---------------------------------------------------------------- design ----------
def test_one_cross_BY_RM_114_segregants(D, ds):
    X = [r for r in rows(D["crosses"]) if r["ds"] == ds]
    assert len(X) == 1 and X[0]["n"] == 114
    P = D["parents"]["name"]
    assert (P[X[0]["p1"]], P[X[0]["p2"]]) == ("BY", "RM")


def test_all_linkages_are_in_that_cross(Q, D, ds):
    x = [i for i, d in enumerate(D["crosses"]["ds"]) if d == ds][0]
    assert all(r["x"] == x and r["m"] == "w" for r in Q)


def test_34_compounds_with_a_linkage(T):
    # "Of the 74 compounds tested, 34 showed at least one significant linkage"
    assert len(T) == 34
    assert all(t["cat"] == "metabolite" for t in T.values())


@pytest.mark.xfail(strict=True, reason=(
    f"Text and the {S1} caption: 'All 52 linkages'. The table has 53 rows (hexose-phosphate and sedoheptulose "
    "7-phosphate each have four)."))
def test_52_linkages(Q):
    assert len(Q) == 52


def test_linkage_count_as_observed(Q):
    assert len(Q) == 53


def test_compounds_with_one_or_two_linkages(Q):
    # "21 of 34 had one mQTL, 9 had two mQTLs"
    c = collections.Counter(per_compound(Q).values())
    assert c[1] == 21 and c[2] == 9


@pytest.mark.xfail(strict=True, reason=(
    f"Text: 'three had three mQTLs and one had four'. {S1}, column Name: two compounds have three rows "
    "(ribose-phosphate, valine) and two have four (hexose-phosphate, Sedoheptulose 7-phosphate)."))
def test_three_with_three_one_with_four(Q):
    c = collections.Counter(per_compound(Q).values())
    assert c[3] == 3 and c[4] == 1


def _parent_sig(Q, raw, T):
    sig = collections.defaultdict(bool)
    for r in Q:
        p = raw[r["id"]][5]
        sig[r["t"]] |= p != "N/A" and float(p) < 0.0898   # caption: FDR 5% corresponds to p = 0.0898
    return sig


@pytest.mark.xfail(strict=True, reason=(
    f"Text: 29 of the 34 compounds with mQTL differed between the parents at FDR 5%. {S1}, column Parent's "
    "P-value, with the caption's cutoff p < 0.0898 on any of a compound's rows: 28 (the rest are N/A or larger)."))
def test_29_of_34_differ_between_parents(Q, raw, T):
    assert sum(_parent_sig(Q, raw, T).values()) == 29


def test_parent_difference_count_as_observed(Q, raw, T):
    assert sum(_parent_sig(Q, raw, T).values()) == 28


def test_lod_above_lowest_threshold(Q):
    # Methods: per-metabolite thresholds range from LOD 3.14 to 3.58
    assert min(r["lod"] for r in Q) >= 3.14


def test_chromosomes_and_cm(Q, D):
    assert all(r["c"] in D["chroms"] for r in Q)
    assert all(r["cm"] is not None and 0 <= r["cm"] < 500 for r in Q)


def test_no_physical_positions(Q):
    """Table S1 has cM only; nothing is placed on R64 (OPEN_QUESTIONS Q1)."""
    assert all(r["a"] is None and r["e"] is None and r["p"] is None for r in Q)


# ---------------------------------------------------------------- hot spots and genes
def test_ira2_hotspot_five_metabolites_on_chrXV(Q, T):
    names = {T[r["t"]]["name"] for r in Q if r["c"] == "chrXV" and r["g"] >= 0}
    assert names == {"citrate/isocitrate", "dihydroxy-acetone phosphate", "Sedoheptulose 7-phosphate",
                     "fructose-1,6-bisphosphate", "hexose-phosphate"}


def test_ura3_hotspot_five_metabolites_on_chrV(Q, T):
    names = {T[r["t"]]["name"] for r in Q if r["c"] == "chrV" and r["g"] >= 0}
    assert names == {"orotate", "orotidine", "orotidine-5'-phosphate", "UDP-D-glucose", "UDP-N-acetyl-glucosamine"}


def test_chrVIII_hotspot_thiamine_SAM_SAH(Q, T):
    names = {T[r["t"]]["name"] for r in Q if r["c"] == "chrVIII"}
    assert names == {"thiamine", "S-adenosyl-L-methionine", "S-adenosyl homocysteine"}


def test_candidate_genes_on_their_linkage_chromosome(D, ds):
    G = D["genes"]
    C = [r for r in rows(D["causal"]) if r["ds"] == ds]
    assert {r["name"] for r in C} == {"IRA2", "URA3", "SLT2", "ERC1"}
    for r in C:
        assert G["c"][r["g"]] == D["qtl"]["c"][r["q1"]] and r["tier"] == "c"


def test_ira2_direction_BY_up(Q, D):
    P = D["parents"]["name"]
    ira = [r for r in Q if r["c"] == "chrXV" and r["g"] >= 0]
    assert ira and all(P[r["up"]] == "BY" for r in ira)


def test_directions_only_where_stated(Q):
    assert sum(1 for r in Q if r["up"] >= 0) == 8


# ---------------------------------------------------------------- provenance ------
def test_every_row_has_resolvable_provenance(D, ds):
    codes = {c for c, d in zip(D["sources"]["code"], D["sources"]["ds"]) if d == ds}
    for t in ("qtl", "causal"):
        for r in rows(D[t]):
            if r["ds"] == ds:
                code, row = r["pv"].split("#")
                assert code in codes and int(row) >= 2


def test_original_rows_kept_verbatim(Q, raw):
    assert len(raw) == len(Q) == 53
    assert all(len(v) == 6 for v in raw.values())
    assert {r["id"] for r in Q} == set(raw)
    for r in Q:
        v = raw[r["id"]]
        assert float(v[2]) == r["lod"] and float(v[4]) == r["cm"]
