"""Validate the ingested Bloom et al. 2019 tables against the paper and against themselves.

    python3 -I -m pytest -q tests/test_bloom2019.py      (after ./build.sh)

Claims from the paper are asserted as stated. Where the data disagree, the test is marked
xfail(strict=True) with the file, column and observed value in the reason: the data are
not changed to match, and if a rebuild ever made the numbers agree the strict xfail would
fail and flag it. Numbers recorded by the authors as comments in the analysis code
(joshsbloom/yeast-16-parents, analysis/jointModelAF.R) are checked the same way.
"""
import collections
import json
import os
import re
import statistics

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "sources", "bloom2019", "raw")
F3D1 = "elife_poa_e49212_Figure_3_source_data_1.xls"
F4D1 = "elife_poa_e49212_Figure_4_source_data_1.xls"
F1D1 = "elife_poa_e49212_Figure_1_source_data_1.xls"
R64 = {"chrI": 230218, "chrII": 813184, "chrIII": 316620, "chrIV": 1531933, "chrV": 576874, "chrVI": 270161,
       "chrVII": 1090940, "chrVIII": 562643, "chrIX": 439888, "chrX": 745751, "chrXI": 666816, "chrXII": 1078177,
       "chrXIII": 924431, "chrXIV": 784333, "chrXV": 1091291, "chrXVI": 948066}
CAMP = ["IRA1", "IRA2", "BCY1", "CYR1", "RAS1", "GPR1"]


def rows(tab):
    keys = list(tab)
    return [dict(zip(keys, v)) for v in zip(*(tab[k] for k in keys))]


@pytest.fixture(scope="module")
def D():
    with open(os.path.join(ROOT, "docs", "data.json"), encoding="utf-8") as f:
        return json.load(f)


@pytest.fixture(scope="module")
def ds(D):
    return D["datasets"]["id"].index("bloom2019")


@pytest.fixture(scope="module")
def Q(D, ds):
    return [r for r in rows(D["qtl"]) if r["ds"] == ds]


@pytest.fixture(scope="module")
def J(Q):
    return [r for r in Q if r["m"] == "j"]


@pytest.fixture(scope="module")
def W(Q):
    return [r for r in Q if r["m"] == "w"]


@pytest.fixture(scope="module")
def CG(D, ds):
    return [r for r in rows(D["causal"]) if r["ds"] == ds]


@pytest.fixture(scope="module")
def genes(D):
    return rows(D["genes"])


# ---------------------------------------------------------------- design ----------
def test_build_is_R64(D):
    assert D["chroms"] == R64
    assert "R64" in D["build"] and "sacCer3" in D["build"]


def test_16_parents(D, ds):
    assert sum(1 for d in D["parents"]["ds"] if d == ds) == 16


def test_16_crosses_each_parent_in_two(D, ds):
    X = [r for r in rows(D["crosses"]) if r["ds"] == ds]
    assert len(X) == 16
    c = collections.Counter([x["p1"] for x in X] + [x["p2"] for x in X])
    assert set(c.values()) == {2}


def test_segregants_13950_avg_872(D, ds):
    n = [r["n"] for r in rows(D["crosses"]) if r["ds"] == ds]
    assert sum(n) == 13950
    assert round(sum(n) / len(n)) == 872


def test_38_traits(D, ds):
    assert sum(1 for d in D["traits"]["ds"] if d == ds) == 38


def test_every_trait_has_a_category(D, ds):
    T = [r for r in rows(D["traits"]) if r["ds"] == ds]
    assert all(t["cat"] and t["catb"] for t in T)


@pytest.mark.skip(reason="298,979 variants: no variant list or genotype matrix in the eLife bundle; not checkable from the files")
def test_298979_variants():
    pass


def test_genotype1_parent_is_alphabetically_first(D, ds):
    """Inferred coding (OPEN_QUESTIONS Q2): consistent with VCF sample order."""
    P = D["parents"]
    for x in rows(D["crosses"]):
        if x["ds"] != ds:
            continue
        assert P["code"][x["g1"]] == min(P["code"][x["p1"]], P["code"][x["p2"]])


# ---------------------------------------------------------------- joint model -----
def test_joint_qtl_4552(J):
    assert len(J) == 4552


def test_joint_per_trait_avg120_range_52_195(J):
    n = collections.Counter(r["t"] for r in J)
    assert len(n) == 38
    assert round(sum(n.values()) / len(n)) == 120
    assert min(n.values()) == 52 and max(n.values()) == 195


def test_joint_per_cross_effects_cover_every_joint_qtl(D, J):
    ids = {i for i, r in enumerate(rows(D["qtl"])) if r["m"] == "j"}
    assert set(D["jeff"]["q"]) <= ids and len(set(D["jeff"]["q"])) == len(J)
    assert len(D["jeff"]["q"]) == 18758


def test_code_comment_common_lead_variants_3568(J):
    # jointModelAF.R: sum(r2MeanEffect$maf1012>.01)  # 3568 / 4552 = 78%
    assert sum(1 for r in J if r["af"] is not None and r["af"] > 0.01) == 3568


def test_code_comment_small_effects_2911(J):
    # jointModelAF.R: sum(abs(r2MeanEffect$MeanAbsBeta)<.1)  # 2911 / 4552 = 64%
    assert sum(1 for r in J if abs(r["ef"]) < 0.1) == 2911


def test_code_comment_effect_by_frequency_table(J):
    # jointModelAF.R:  rare <.3 839, rare >.3 145, common <.3 3478, common >.3 90
    t = collections.Counter((r["af"] > 0.01, abs(r["ef"]) > 0.3) for r in J)
    assert (t[(False, False)], t[(False, True)], t[(True, False)], t[(True, True)]) == (839, 145, 3478, 90)


def test_code_comment_median_abs_effect(J):
    # jointModelAF.R: median(abs(r2MeanEffect$MeanAbsBeta))  # 0.08071
    assert round(statistics.median(abs(r["ef"]) for r in J), 5) == pytest.approx(0.08071, abs=1e-5)


# ---------------------------------------------------------------- within-cross ----
def test_within_cross_all_q_below_005(W):
    assert len(W) == 7751
    assert max(r["q"] for r in W) < 0.05


def test_within_cross_median_12_per_trait(W, D):
    per = collections.Counter((r["x"], r["t"]) for r in W)
    assert statistics.median(per.values()) == 12
    a = D["crosses"]["code"].index("A")
    # mapping.R reports the BYxRM cross (A) separately; the other 15 crosses give 12 too
    assert statistics.median(v for (x, t), v in per.items() if x != a) == 12


def _joint_only(J, W, same_cross, D):
    jeff = collections.defaultdict(set)
    for q, x in zip(D["jeff"]["q"], D["jeff"]["x"]):
        jeff[D["qtl"]["id"][q]].add(x)
    iv = collections.defaultdict(list)
    for r in W:
        iv[(r["t"], r["c"])].append((r["a"], r["e"], r["x"]))
    n = 0
    for r in J:
        hits = [x for a, e, x in iv[(r["t"], r["c"])] if a <= r["p"] <= e]
        if same_cross:
            hits = [x for x in hits if x in jeff[r["id"]]]
        n += not hits
    return n


@pytest.mark.xfail(strict=True, reason=(
    "Paper: 458 joint QTL not detected within crosses. Data: 1,022 joint QTL whose lead variant falls in no "
    "1.5-LOD interval of the same trait in a cross where it segregates (the jointModelAF.R criterion), 526 if any "
    "cross counts. Sources: Figure 3-source data 1, sheets 'joint model (effect per cross)' (peaks, cross) and "
    "'within-cross model' (1.5 LOD drop CI left/right). The code overlapped against cross.peaks, which holds "
    "within-cross QTL up to q < 0.2; the published sheet holds only q < 0.05, so the count is not reproducible "
    "from the bundle."))
def test_joint_only_458(J, W, D):
    assert _joint_only(J, W, True, D) == 458


def test_joint_only_counts_as_observed(J, W, D):
    assert _joint_only(J, W, True, D) == 1022
    assert _joint_only(J, W, False, D) == 526


# ---------------------------------------------------------------- coordinates -----
def test_intervals_inside_chromosomes(Q):
    for r in Q:
        assert r["c"] in R64
        assert 1 <= r["a"] and r["e"] <= R64[r["c"]], r["id"]


def test_start_le_peak_le_end(Q):
    bad = [r["id"] for r in Q if not (r["a"] <= r["p"] <= r["e"])]
    assert not bad


def test_joint_intervals_are_the_lead_variant(J):
    assert all(r["a"] == r["p"] == r["e"] and r["lv"].split("_")[0] == r["c"] and int(r["lv"].split("_")[1]) == r["p"] for r in J)


# ---------------------------------------------------------------- candidate genes -
def test_cAMP_genes_in_causal_table(CG):
    names = {r["name"] for r in CG}
    assert set(CAMP) <= names


@pytest.mark.parametrize("gene", CAMP)
def test_cAMP_gene_at_fdr5(CG, gene):
    if gene == "BCY1":
        pytest.xfail("BCY1's only Causal_Genes row (Figure 4-source data 1, column localFDR) has local FDR 0.1507: "
                     "FDR 20% tier, not FDR 5%")
    assert min(r["fdr"] for r in CG if r["name"] == gene) < 0.05


@pytest.mark.xfail(strict=True, reason=(
    "Paper: 105 unique candidate causal genes at FDR 5%. Data: 102 unique ORFs with localFDR < 0.05 "
    "(also 102 with <= 0.05) in Figure 4-source data 1, sheet Causal_Genes, columns ORF and localFDR; "
    "195 unique ORFs in the whole sheet (local FDR < 0.2)."))
def test_105_genes_at_fdr5(CG):
    assert len({r["orf"] for r in CG if r["fdr"] < 0.05}) == 105


def test_gene_counts_as_observed(CG):
    assert len({r["orf"] for r in CG if r["fdr"] < 0.05}) == 102
    assert len({r["orf"] for r in CG}) == 195


def test_previously_identified_genes_all_in_causal_table(D, ds):
    names = {r["name"] for r in rows(D["causal"]) if r["ds"] == ds}
    with open(os.path.join(ROOT, "sources", "bloom2019", "code_tables", "previously_identified.tsv")) as f:
        prev = [l.split("\t")[0] for l in f if not l.startswith("#")][1:]
    assert len(prev) == 25 and set(prev) <= names


@pytest.mark.xfail(strict=True, reason=(
    "Paper: 24 of 26 experimentally validated genes found at FDR 5%. The bundle has no list of the 26. "
    "The analysis code's list (QTL_causality.R lines 438-454; sources/bloom2019/code_tables/previously_identified.tsv) "
    "has 25 names; 23 have a Causal_Genes row with localFDR < 0.05 and HAP1, RGA1 only < 0.2."))
def test_24_validated_at_fdr5(CG):
    with open(os.path.join(ROOT, "sources", "bloom2019", "code_tables", "previously_identified.tsv")) as f:
        prev = {l.split("\t")[0] for l in f if not l.startswith("#")} - {"name"}
    assert len(prev & {r["name"] for r in CG if r["fdr"] < 0.05}) == 24


def _overlap(a, e, lo, hi):
    return a <= hi and e >= lo


def test_assigned_gene_within_code_search_window(D, W, CG):
    """Every within-cross QTL's assigned gene lies in the window the code fine-mapped:
    25 kb either side of the peak (QTL_causality.R), using the sheet's gene region."""
    by_q = collections.defaultdict(list)
    for r in CG:
        for q in (r["q1"], r["q2"]):
            if q >= 0:
                by_q[q].append(r)
    ids = {r["id"]: i for i, r in enumerate(rows(D["qtl"]))}
    n = 0
    for r in W:
        if r["g"] < 0:
            continue
        rs = [c for c in by_q[ids[r["id"]]] if c["g"] == r["g"]]
        assert rs, r["id"]
        assert any(c["c"] == r["c"] and _overlap(c["a"], c["e"], r["p"] - 25000, r["p"] + 25000) for c in rs), r["id"]
        n += 1
    assert n > 0


def _within_outside(W, genes):
    out = []
    for r in W:
        if r["g"] < 0:
            continue
        g = genes[r["g"]]
        if not (g["c"] == r["c"] and _overlap(g["a"], g["e"], r["a"], r["e"])):
            out.append(r["id"])
    return out


@pytest.mark.xfail(strict=True, reason=(
    "Some candidate genes lie outside the 1.5-LOD interval of the within-cross QTL they are assigned to "
    "(SGD R64-5-1 gene bounds vs Figure 3-source data 1 'within-cross model' CI columns; the gene comes from "
    "Figure 4-source data 1 Causal_Genes). The code searched 25 kb either side of each peak, which can extend "
    "past the 1.5-LOD interval. Counts are printed by test_assigned_gene_interval_counts."))
def test_assigned_gene_inside_lod_interval(W, genes):
    assert not _within_outside(W, genes)


def test_assigned_gene_interval_counts(W, genes, capsys):
    assigned = [r for r in W if r["g"] >= 0]
    out = _within_outside(W, genes)
    with capsys.disabled():
        print(f"\n  within-cross QTL with a candidate gene: {len(assigned)}; gene outside the 1.5-LOD interval: {len(out)}")
    assert 0 < len(out) < len(assigned)


def test_joint_assigned_gene_near_lead_variant(J, genes, capsys):
    """Joint QTL: the code marks a lead variant resolved when it is one of the fine-mapped
    variants of a resolved QTL pair (jointModelAF.R). Both the variant and the gene lie
    within 25 kb of that pair's peaks, so they are at most 50 kb apart."""
    res = [r for r in J if r["g"] >= 0]
    assert len(res) == 385
    inside = 0
    for r in res:
        g = genes[r["g"]]
        assert g["c"] == r["c"]
        d = max(0, g["a"] - r["p"], r["p"] - g["e"])
        assert d <= 50000, r["id"]
        inside += d == 0
    with capsys.disabled():
        print(f"\n  joint QTL with a candidate gene: {len(res)}; lead variant inside the SGD gene: {inside}")


def test_tiers(D, ds):
    C = [r for r in rows(D["cand"]) if r["ds"] == ds]
    assert {c["tier"] for c in C} == {"v", "5", "20"}
    assert all((c["tier"] == "5") == (c["fdr"] < 0.05) for c in C if c["tier"] != "v")


# ---------------------------------------------------------------- source cross-checks
def test_supp_file_3_within_cross_counts(W, D):
    hdr = _supp3()
    T = D["traits"]
    per = collections.Counter(T["raw"][r["t"]] for r in W)
    for raw, (total, _joint) in hdr.items():
        assert per[raw] == total, raw


@pytest.mark.xfail(strict=True, reason=(
    "Supplementary file 3 page headers ('N total QTL | M joint QTL') give one more joint QTL than "
    "Figure 3-source data 1, sheet 'joint model (average effect)', for EtOH_Glucose;;1 (111 vs 110), "
    "Glycerol;;1 (116 vs 115) and Lithium_Chloride;100mM;2 (159 vs 158). All other 35 traits agree."))
def test_supp_file_3_joint_counts(J, D):
    hdr = _supp3()
    per = collections.Counter(D["traits"]["raw"][r["t"]] for r in J)
    assert {raw: m for raw, (_, m) in hdr.items()} == {raw: per[raw] for raw in hdr}


def test_supp_file_3_joint_mismatches_are_exactly_three(J, D):
    hdr = _supp3()
    per = collections.Counter(D["traits"]["raw"][r["t"]] for r in J)
    bad = {raw: (m, per[raw]) for raw, (_, m) in hdr.items() if per[raw] != m}
    assert bad == {"EtOH_Glucose;;1": (111, 110), "Glycerol;;1": (116, 115), "Lithium_Chloride;100mM;2": (159, 158)}


_SUPP3 = {"6−azauracil": "6-azauracil;50ug/mL;3", "EtOH Glu": "EtOH_Glucose;;1", "YNB ph3": "YNB;ph3;1",
          "YNB ph8": "YNB;ph8;1", "YPD": "YPD;;1", "YPD 15C": "YPD;15;1", "YPD 37C": "YPD;37;1", "YNB": "YNB;;1"}


def _supp3():
    pdfplumber = pytest.importorskip("pdfplumber")
    out = {}
    with pdfplumber.open(os.path.join(RAW, "elife_poa_e49212_Supplementary_file_3.pdf")) as pdf:
        for pg in pdf.pages:
            m = re.search(r"([^\n]*?) (\d+) total QTL \| (\d+) joint QTL", pg.extract_text())
            name = m.group(1).strip()
            raw = _SUPP3.get(name)
            if raw is None:
                comp = name.replace(" ", "_")
                raw = {"Copper_Sulfate": "Copper_Sulfate;;1"}.get(comp)
            out[raw or name] = (int(m.group(2)), int(m.group(3)))
    # map the remaining display names onto trait strings by compound
    return _resolve(out)


def _resolve(out):
    import xlrd
    s = xlrd.open_workbook(os.path.join(RAW, F3D1)).sheet_by_name("joint model (average effect)")
    traits = {s.cell_value(i, 0) for i in range(1, s.nrows)}
    res = {}
    for k, v in out.items():
        if k in traits:
            res[k] = v
            continue
        hit = [t for t in traits if t.split(";")[0].replace("_", " ").lower() == k.lower()]
        assert len(hit) == 1, k
        res[hit[0]] = v
    assert len(res) == 38
    return res


# ---------------------------------------------------------------- provenance ------
def test_every_row_has_resolvable_provenance(D):
    codes = set(D["sources"]["code"])
    for t in ("crosses", "qtl", "jeff", "causal", "vexp"):
        for pv in D[t]["pv"]:
            code, row = pv.split("#")
            assert code in codes and int(row) >= 2


def test_original_rows_kept_verbatim(D):
    """Every QTL, per-cross effect and causal row has its full source row in docs/v/."""
    sc = D["datasets"]["sc"][D["datasets"]["id"].index("bloom2019")]
    seen = collections.Counter()
    width = {}
    for ch in R64:
        with open(os.path.join(ROOT, "docs", "v", f"{sc}_{ch}.json"), encoding="utf-8") as f:
            v = json.load(f)
        for code, rs in v["rows"].items():
            width[code] = len(v["cols"][code])
            for rid, val in rs.items():
                if code == "F3D1:pc":
                    seen[code] += len(val)
                    assert all(len(x) == width[code] for x in val)
                else:
                    seen[code] += 1
                    assert len(val) == width[code]
    assert seen == {"F3D1:j": 4552, "F3D1:w": 7751, "F3D1:pc": 18758, "F4D1:cg": 427}
    assert width == {"F3D1:j": 11, "F3D1:w": 14, "F3D1:pc": 9, "F4D1:cg": 66}


def test_variance_explained_table(D):
    assert len(D["vexp"]["t"]) == 608 == 16 * 38
