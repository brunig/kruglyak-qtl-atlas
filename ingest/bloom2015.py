"""Ingest Bloom et al. 2015 (Nature Communications 6:8712), 4,390 BY x RM segregants, 20 traits.

    python3 -I ingest/bloom2015.py          # from the repository root

Reads only from sources/: the supplementary files in sources/bloom2015/raw/ (Supplementary Data
1-4 as .xls), sources/bloom2015/extracted/markers.tsv (the genotype-matrix column names of
Supplementary Data 5, which the `index` columns of Supplementary Data 3 and 4 refer to, 1-based)
and supp_data_5_contents.tsv, and the SGD tables in sources/sgd/. Writes
build/fragments/bloom2015.json, docs/v/b15_*.json and, through ingest/atlas.py, docs/data.json.
Deterministic. No genotype or phenotype value from Supplementary Data 5 is used or copied.

Genome build: S288C R64. The Methods align reads to sacCer3; every one of the 28,220 marker
positions lies inside its R64 chromosome, and every Supplementary Data 3 `pos` equals the
position of the marker its `index` names. No coordinate transform.

Semantics (article Methods; OPEN_QUESTIONS.md, Bloom 2015 section):
- Additive QTL: forward stepwise scan per chromosome on BLUP residuals, 5% FWER from 1,000
  permutations; `pos` is the forward-search peak marker. Intervals are 1.5-LOD drops from R/qtl
  (`pos.left/peak/right_refineqtl`); `max.lod` is read as the maximum LOD of that refineqtl
  profile (most defensible reading; some values are below 3).
- lm_coefficient: coefficient of the peak marker in a per-trait multiple regression of the
  phenotype on all significant peak markers, genotypes coded -1 (BY) / 1 (RM) as Supplementary
  Data 5 states, so a positive value means the RM allele increases growth (most defensible
  reading). Phenotypes are normalized end-point colony radius, centred near 0 but not scaled to
  unit variance, so the coefficient is in phenotype units per genotype step (half the allele
  difference).
- QTL whose interval holds exactly one gene are the paper's single-gene QTL; that gene is
  recorded as a candidate (tier c).
"""
import collections
import hashlib
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import atlas  # noqa: E402

DS_ID, SC = "bloom2015", "b15"
SRC = os.path.join(atlas.SOURCES, "bloom2015")
RAW = os.path.join(SRC, "raw")
P = "41467_2015_BFncomms9712_MOESM7"
F = {"SI": P + "24_ESM.pdf", "SD1": P + "25_ESM.xls", "SD2": P + "26_ESM.xls", "SD3": P + "27_ESM.xls",
     "SD4": P + "28_ESM.xls", "SD5": P + "29_ESM.zip"}
MK = "../extracted/markers.tsv"                # paths relative to raw/, as sources.file is
SD5L = "../extracted/supp_data_5_contents.tsv"
ART = "../article.html"

# The paper does not group its traits; categories follow the same labels as the other datasets.
CHEM = "atlas assignment: compound added to the medium (dose not given in this paper; see OPEN_QUESTIONS)"
CARB = "atlas assignment: alternative carbon source"
CATEGORY = {
    "CobaltChloride": ("drug/chemical", CHEM), "CopperSulfate": ("drug/chemical", CHEM), "Diamide": ("drug/chemical", CHEM),
    "E6-Berbamine": ("drug/chemical", CHEM), "Formamide": ("drug/chemical", CHEM), "Hydroxyurea": ("drug/chemical", CHEM),
    "IndolaceticAcid": ("drug/chemical", CHEM), "MagnesiumChloride": ("drug/chemical", CHEM),
    "ManganeseSulfate": ("drug/chemical", CHEM), "Menadione": ("drug/chemical", CHEM), "Neomycin": ("drug/chemical", CHEM),
    "Zeocin": ("drug/chemical", CHEM),
    "Lactate": ("carbon source", CARB), "Lactose": ("carbon source", CARB), "Raffinose": ("carbon source", CARB),
    "Trehalose": ("carbon source", CARB), "Xylose": ("carbon source", CARB),
    "Ethanol": ("carbon source", "atlas assignment, ambiguous: the paper does not say whether ethanol replaced glucose or was added to it"),
    "YNB": ("reference medium", "atlas assignment: minimal medium without additions"),
    "YPD": ("reference medium", "atlas assignment: rich medium without additions"),
}


def num(v):
    if isinstance(v, float) and v.is_integer() and abs(v) < 1e15:
        return int(v)
    return v


def verbatim(row):
    return [None if v == "" else num(v) for v in row]


def sheet(name, idx_or_name):
    import xlrd
    b = xlrd.open_workbook(os.path.join(RAW, F[name]))
    s = b.sheet_by_name(idx_or_name) if isinstance(idx_or_name, str) else b.sheet_by_index(idx_or_name)
    return s.name, [s.row_values(i) for i in range(s.nrows)]


def manifest():
    out = {}
    with open(os.path.join(SRC, "MANIFEST.tsv")) as f:
        for line in f.read().splitlines()[1:]:
            name, _, sha = line.split("\t")
            out[name] = sha
    for rel in (MK, SD5L, ART):
        out[rel] = hashlib.sha256(open(os.path.normpath(os.path.join(RAW, rel)), "rb").read()).hexdigest()
    return out


def markers():
    out = {}
    with open(os.path.join(SRC, "extracted", "markers.tsv")) as f:
        rows = [l.rstrip("\n").split("\t") for l in f if not l.startswith("#")]
    for r in rows[1:]:
        out[int(r[0])] = (r[2], int(r[3]), r[4], r[5])
    return out


def display(raw):
    if raw in ("YNB", "YPD"):
        return raw
    if raw == "E6-Berbamine":
        return "E6-berbamine"
    words = re.findall(r"[A-Z][a-z0-9]*|[a-z0-9]+", raw)
    s = " ".join(words).lower()
    return s[0].upper() + s[1:]


def slug(s):
    return re.sub(r"-+", "-", re.sub(r"[^A-Za-z0-9]+", "-", s)).strip("-")


def main():
    sha = manifest()
    MKR = markers()
    C = collections.OrderedDict
    T = {k: C() for k in atlas.TABLES}

    def add(table, **kw):
        tab = T[table]
        n = atlas.nlen(tab)
        for k, v in kw.items():
            if k not in tab:
                tab[k] = [None] * n
            tab[k].append(v)
        for k in tab:
            if k not in kw:
                tab[k].append(None)
        return n

    s3name, s3 = sheet("SD3", 0)
    f4name, s4f = sheet("SD4", "Full Scan")
    m4name, s4m = sheet("SD4", "Marginal Scan")
    s2name, s2 = sheet("SD2", 0)
    s1name, s1 = sheet("SD1", 0)
    srcs = {"b15:SD3": (F["SD3"], s3name, 2), "b15:SD4f": (F["SD4"], f4name, 1), "b15:SD4m": (F["SD4"], m4name, 1),
            "b15:SD2": (F["SD2"], s2name, 3), "b15:SD1": (F["SD1"], s1name, 3), "b15:SD5": (F["SD5"], "(zip archive)", None),
            "b15:MK": (MK, "markers", 3), "b15:SD5L": (SD5L, "contents", 3), "b15:ART": (ART, "article text", None)}
    for code, (f, sh, hrow) in srcs.items():
        add("sources", code=code, ds=0, file=f, sheet=sh, hrow=hrow, sha256=sha[f])

    # ---- parents and the one cross (strain descriptions from the article's Introduction) ---
    by = add("parents", id=f"{SC}-BY", ds=0, name="BY", code="BY", s1=None, s2=None, p18=None, org="laboratory strain (article text)")
    rm = add("parents", id=f"{SC}-RM", ds=0, name="RM", code="RM", s1=None, s2=None, p18=None, org="vineyard isolate (article text)")
    with open(os.path.join(SRC, "extracted", "supp_data_5_contents.tsv")) as f:
        geno_line = [i + 1 for i, l in enumerate(f) if "\tgdata\t" in l][0]
    # genotypes are coded -1 (BY) / 1 (RM): BY is the lower code, recorded as genotype 1 (crosses.g1)
    x = add("crosses", id=f"{SC}-BYxRM", ds=0, code="BYxRM", name="BY × RM", p1=by, p2=rm, g1=by, n=4390,
            dip=None, mm=None, pv=f"b15:SD5L#{geno_line}")

    # ---- traits ------------------------------------------------------------------------------
    head = [str(v) for v in s3[1]]
    assert head[:6] == ["chr", "trait", "index", "pos", "lm_coefficient", "variance_explained"], head
    assert head[7:14] == ["pos.left_refineqtl", "pos.peak_refineqtl", "pos.right_refineqtl", "max.lod", "int.size",
                          "systematic name", "common name"], head
    data = [(i + 1, r) for i, r in enumerate(s3) if i >= 2 and r[0] != ""]
    tidx = {}
    for raw in sorted({r[1] for _, r in data}, key=str.lower):
        cat, catb = CATEGORY[raw]
        tidx[raw] = add("traits", id=f"{SC}-{slug(raw)}", ds=0, raw=raw, name=display(raw), cond=None, bt=None,
                        ck=display(raw).lower(), cat=cat, catb=catb)

    G = atlas.genes()
    gidx = {}
    for i, g in enumerate(G["id"]):
        gidx.setdefault(g, i)

    orig = collections.defaultdict(lambda: collections.defaultdict(dict))
    for rn, r in data:
        chrom, trait, idx, pos = r[0], r[1], int(r[2]), int(r[3])
        mc, mp, mref, malt = MKR[idx]
        assert (mc, mp) == (chrom, pos), (rn, chrom, pos, MKR[idx])
        ef = float(r[4])
        genes = str(r[12]).split()
        g = gidx.get(genes[0], -1) if len(genes) == 1 else -1
        qid = f"{SC}-w-{rn}"
        q = add("qtl", id=qid, ds=0, t=tidx[trait], x=x, m="w", c=chrom, a=int(r[7]), e=int(r[9]), p=int(r[8]),
                lod=atlas.r6(r[10]), q=None, ef=atlas.r6(ef), ve=atlas.r6(r[5]), lv=f"{mc}_{mp}_{mref}_{malt}",
                af=None, cc=None, pa=None, anc=None, daf=None, up=rm if ef > 0 else by, ups=None, dns=None, xs=None,
                k=None, g=g, gf=None, gt="c" if g >= 0 else None, pv=f"b15:SD3#{rn}")
        orig[chrom]["b15:SD3"][qid] = verbatim(r)
        if g >= 0:
            cid = f"{SC}-c-{rn}"
            add("causal", id=cid, ds=0, t=tidx[trait], g=g, orf=G["id"][g], name=str(r[13]).strip() or None,
                c=G["c"][g], a=G["a"][g], e=G["e"][g], pcs=None, sp=None, eff=None, fdr=None, jppc=None, jv=None,
                q1=q, q2=-1, tier="c", pv=f"b15:SD3#{rn}")
            orig[chrom]["b15:SD3"][cid] = verbatim(r)
    cz = T["causal"]
    per_gene = C()
    for i in range(atlas.nlen(cz)):
        per_gene.setdefault(cz["g"][i], []).append(i)
    for g, ii in per_gene.items():
        add("cand", ds=0, g=g, tier="c", fdr=None, nt=len({cz["t"][i] for i in ii}), n=len(ii))

    nq = atlas.nlen(T["qtl"])
    add("datasets", id=DS_ID, sc=SC, label="Bloom 2015",
        cite="Bloom JS, Kotenko I, Sadhu MJ, Treusch S, Albert FW, Kruglyak L. 2015. Genetic interactions contribute less than additive effects to quantitative trait variation in yeast. Nature Communications 6:8712",
        doi="10.1038/ncomms9712", year=2015,
        design="segregant panel; 4,390 BY × RM segregants (an expansion of the Bloom 2013 panel), 20 end-point growth traits with two replicates; additive QTL and QTL-QTL interactions",
        build="S288C R64 (reads aligned to sacCer3, Methods). All 28,220 marker positions lie inside R64 chromosomes and each QTL's pos equals the position of the marker its index names. No transform.",
        method="Additive QTL: forward stepwise scan per chromosome on BLUP residuals (mixed model with the rest of the genome and pairwise interactions as random effects), LOD from the marker-residual correlation, 5% FWER from 1,000 permutations, repeated until no further QTL. "
               "QTL-QTL interactions: full two-dimensional and marginal (additive QTL × genome) scans on BLUP residuals, FDR 10% from five permutations; kept in docs/v/b15_extra.json, not drawn.",
        fdr=None,
        interval="1.5-LOD drop interval from R/qtl (pos.left/right_refineqtl); qtl.p is the refineqtl peak (pos.peak_refineqtl). qtl.lv is the forward-search peak marker (Supplementary Data 3 pos/index), which can differ from p and in 11 QTL lies outside the interval. qtl.lod is max.lod of the refineqtl profile.",
        effect="lm_coefficient: coefficient of the peak marker in a per-trait regression on all significant peak markers, genotypes coded -1 (BY) / 1 (RM); positive = the RM allele increases growth (most defensible reading). Units: normalized colony radius (centred near 0, not unit variance) per genotype step, i.e. half the allele difference.",
        nq=nq,
        notes="1 cross, 4,390 segregants, 20 traits. Additive QTL at 5% FWER per chromosome (Methods), so datasets.fdr is empty. "
              "Candidate genes: the QTL whose interval contains a single gene (tier c), as listed in Supplementary Data 3; not claimed causal beyond the paper. "
              "205 full-scan and 266 marginal-scan QTL-QTL interactions are kept verbatim, not drawn.")

    atlas.save_fragment(DS_ID, T)

    # ---- lazy files ------------------------------------------------------------------------------
    s3cols = [str(v) if v != "" else "" for v in s3[1]]
    for chrom in atlas.chrom_lengths():
        out = {"cols": {"b15:SD3": s3cols}, "rows": {"b15:SD3": dict(orig[chrom].get("b15:SD3", {}))}}
        atlas.write_json(os.path.join(atlas.DOCS, "v", f"{SC}_{chrom}.json"), out)

    def keep(rows, hrow):
        return {"cols": [str(v) for v in rows[hrow - 1]], "rows": [verbatim(r) for r in rows]}

    with open(os.path.join(SRC, "extracted", "supp_data_5_contents.tsv")) as f:
        sd5 = [l.rstrip("\n").split("\t") for l in f if not l.startswith("#")]
    extra = {"b15:SD4f": keep(s4f, 1), "b15:SD4m": keep(s4m, 1), "b15:SD2": keep(s2, 3), "b15:SD1": keep(s1, 3),
             "b15:SD5L": {"cols": sd5[0], "rows": sd5[1:]}}
    atlas.write_json(os.path.join(atlas.DOCS, "v", f"{SC}_extra.json"), extra)

    data_ = atlas.assemble()
    print(f"{DS_ID}: 1 cross (4,390 segregants), {len(tidx)} traits, {nq} additive QTL, "
          f"{atlas.nlen(T['causal'])} single-gene QTL ({atlas.nlen(T['cand'])} genes), "
          f"{len(s4f) - 1} full-scan and {len(s4m) - 1} marginal-scan interactions kept verbatim; "
          f"docs/data.json: {len(data_['datasets']['id'])} dataset(s)", file=sys.stderr)


if __name__ == "__main__":
    main()
