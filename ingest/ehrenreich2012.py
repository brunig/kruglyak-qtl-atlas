"""Ingest Ehrenreich et al. 2012 (PLOS Genetics 8:e1002570): X-QTL for 13 chemicals in all six
pairwise crosses of BY, RM, YJM789 and YPS163.

    python3 -I ingest/ehrenreich2012.py          # from the repository root

Reads only from sources/ (sources/ehrenreich2012/raw/, article.html, sources/sgd/). Writes
build/fragments/ehrenreich2012.json, docs/v/e12_*.json and, through ingest/atlas.py,
docs/data.json. Deterministic.

Semantics (File S10 header, File S8/S9 peak caller, article Methods):
- One File S10 row = one distinct locus for one trait, after grouping the peaks of all six crosses
  within a 200 kb window around the most strongly selected peak. `position` is the locus position
  (most defensible reading: the position of that strongest peak; the file does not say).
- `<cross>` TRUE/FALSE: a peak was detected in that cross at 1% FDR. `<cross>_sign`: +1 the
  allele of the first-named (MATalpha) parent was enriched in the selected pool, -1 the
  MATa parent's, 0 no peak. Enrichment in a resistance selection means that parent's allele
  confers resistance. (Confirmed: the BY-resistant singleton pattern 1_1_1_0_0_0 has count 10 in
  both File S10 and Table S1, and Copper sulfate BYxRM has 18 BY / 2 RM loci, as in the text.)
- No interval, LOD or effect size exists; loci are points, direction is by strain.

Coordinates: positions agree with R64 gene positions (CAN1, ENT3/VPS70/RSF2 and MKT1 clusters;
every File S11 gene on the right chromosome within 15 kb), and no position exceeds an R64
chromosome. The SGD release used for array design is not stated; positions are used as given,
with no transform (OPEN_QUESTIONS.md, Ehrenreich 2012 section, Q1).
"""
import ast
import collections
import csv
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import atlas  # noqa: E402

DS_ID, SC = "ehrenreich2012", "e12"
SRC = os.path.join(atlas.SOURCES, DS_ID)
RAW = os.path.join(SRC, "raw")
PFX = "journal.pgen.1002570."
ROMAN = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI"]
CROSSES = ["BYxRM", "BYxYJM", "BYxYPS", "RMxYJM", "RMxYPS", "YJMxYPS"]   # File S10 column order
# Parents, with the descriptions and origins the article gives (Introduction, Author Summary, Methods).
PARENTS = [
    ("BY", "BY", "BY4716 (laboratory strain, S288C background)", "laboratory strain"),
    ("RM", "RM", "RM11-1a (vineyard strain)", "vineyard"),
    ("YJM789", "YJM", "YJM789, derived from a clinical isolate (immunocompromised human)", "clinical isolate"),
    ("YPS163", "YPS", "YPS163, oak strain (sap of an oak tree); arrays designed from the isogenic YPS606 genome", "oak tree sap"),
]
DISPLAY = {"LY": "LY83583", "4-NQO": "4-NQO"}

SOURCES = {   # provenance code -> (file, sheet, 1-based header row)
    "e12:S10": (PFX + "s013.txt", "(text table)", 1),
    "e12:S11": (PFX + "s014.txt", "(text table)", 1),
    "e12:S1dose": (PFX + "s004.xlsx", "Drug doses used in X-QTL", 2),
    "e12:S1det": (PFX + "s004.xlsx", "Determination of selective dose", 1),
}


def manifest():
    out = {}
    with open(os.path.join(SRC, "MANIFEST.tsv")) as f:
        for line in f.read().splitlines()[1:]:
            name, _, sha = line.split("\t")
            out[name] = sha
    return out


def read_tsv(code):
    """Return (column names, [(file line number, values)])."""
    path = os.path.join(RAW, SOURCES[code][0])
    with open(path, encoding="utf-8", newline="") as f:
        rows = list(csv.reader(f, delimiter="\t"))
    return rows[0], [(i + 2, r) for i, r in enumerate(rows[1:])]


def read_xlsx(code):
    import openpyxl
    import warnings
    fname, sheet, hrow = SOURCES[code]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        b = openpyxl.load_workbook(os.path.join(RAW, fname), read_only=True, data_only=True)
    rows = [["" if v is None else v for v in r] for r in b[sheet].iter_rows(values_only=True)]
    width = max(len(r) for r in rows)
    rows = [list(r) + [""] * (width - len(r)) for r in rows]
    return [str(v) for v in rows[hrow - 1]], [(i + 1, r) for i, r in enumerate(rows) if i >= hrow]


def verbatim(row):
    out = []
    for v in row:
        if v in ("", None):
            out.append(None)
        elif isinstance(v, float) and v.is_integer():
            out.append(int(v))
        else:
            out.append(v)
    return out


def slug(s):
    return re.sub(r"-+", "-", re.sub(r"[^A-Za-z0-9]+", "-", s)).strip("-")


def main():
    sha = manifest()
    T = {k: collections.OrderedDict() for k in atlas.TABLES}

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

    for code, (fname, sheet, hrow) in SOURCES.items():
        add("sources", code=code, ds=0, file=fname, sheet=sheet, hrow=hrow, sha256=sha[fname])

    # ---- parents and crosses ------------------------------------------------------
    pidx = {}
    for name, code, desc, org in PARENTS:
        pidx[code] = add("parents", id=f"{SC}-{name}", ds=0, name=name, code=code, s1=desc, s2=desc, p18=None, org=org)
    dcols, drows = read_xlsx("e12:S1dose")
    assert dcols[1:7] == CROSSES, dcols
    xidx = {}
    for c in CROSSES:
        a, b = c.split("x")
        xidx[c] = add("crosses", id=f"{SC}-{c}", ds=0, code=c, name=c, p1=pidx[a], p2=pidx[b], g1=None,
                      n=None, dip=None, mm=None, pv="e12:S1dose#2")

    # ---- traits (13 chemicals; dose is cross-specific, File S1) ---------------------
    dose = {}
    for rn, r in drows:
        if r[0]:
            dose[r[0]] = {c: r[1 + k] for k, c in enumerate(CROSSES)}
    qcols, qrows = read_tsv("e12:S10")
    traits = sorted({r[0] for _, r in qrows}, key=str.lower)
    dose_name = {"LY": "LY83583", "Sodium dodecyl sulfate": "SDS"}   # File S1 spells these two differently
    tidx = {}
    for raw in traits:
        dd = dose[dose_name.get(raw, raw)]
        cond = "dose by cross: " + "; ".join(f"{c} {dd[c]}" for c in CROSSES)
        name = DISPLAY.get(raw, raw[0].upper() + raw[1:])
        ck = {"Sodium dodecyl sulfate": "sds"}.get(raw, DISPLAY.get(raw, raw).lower())
        tidx[raw] = add("traits", id=f"{SC}-{slug(raw)}", ds=0, raw=raw, name=name, cond=cond, bt=None, ck=ck,
                        cat="drug/chemical", catb="File S1 (X-QTL selection doses), drug selection")

    # ---- loci (File S10) --------------------------------------------------------------
    assert qcols == ["trait", "chromosome", "position", "grouping"] + CROSSES + [c + "_sign" for c in CROSSES] + ["pattern"], qcols
    orig = collections.defaultdict(lambda: collections.defaultdict(dict))
    qkey = {}
    for rn, r in qrows:
        trait, chrom, pos = r[0], "chr" + ROMAN[int(r[1]) - 1], int(float(r[2]))
        xs, ups, dns = [], [], []
        for k, c in enumerate(CROSSES):
            det, sgn = r[4 + k], int(r[10 + k])
            assert (det == "TRUE") == (sgn != 0), (rn, c)
            if sgn:
                a, b = c.split("x")
                xs.append(xidx[c]); ups.append(pidx[a if sgn > 0 else b]); dns.append(pidx[b if sgn > 0 else a])
        qid = f"{SC}-x-{rn}"
        q = add("qtl", id=qid, ds=0, t=tidx[trait], x=-1, xs=xs, k="grouped X-QTL locus", m="x", c=chrom, a=pos, e=pos, p=pos,
                lod=None, q=None, ef=None, ve=None, lv=None, af=None, cc=None, pa=None, anc=None, daf=None,
                up=-1, ups=ups, dns=dns, g=-1, gf=None, gt=None, pv=f"e12:S10#{rn}")
        qkey[(trait, r[1], r[2])] = q
        orig[chrom]["e12:S10"][qid] = verbatim(r)

    # ---- candidate genes (File S11) -------------------------------------------------------
    G = atlas.genes()
    gidx = {}
    for i, g in enumerate(G["id"]):
        gidx.setdefault(g, i)
    ccols, crows = read_tsv("e12:S11")
    assert ccols == ["PeakNum", "Trait", "DifferentStrain", "Chromosome", "Position", "WindowSize", "Genes"], ccols
    missing_genes = []
    for rn, r in crows:
        q = qkey[(r[1], r[3], r[4])]
        sp = pidx[r[2]]
        # direction of the singleton strain's allele, from that locus's File S10 signs
        fav = [u == sp for u, d in zip(T["qtl"]["ups"][q], T["qtl"]["dns"][q]) if sp in (u, d)]
        eff = "increase" if fav and all(fav) else "decrease" if fav and not any(fav) else None
        genes = ast.literal_eval(r[6])
        assert len(genes) == int(r[5]), rn
        chrom = "chr" + ROMAN[int(r[3]) - 1]
        for j, orf in enumerate(genes):
            g = gidx.get(orf, -1)
            if g < 0:
                missing_genes.append(orf)
            cid = f"{SC}-c-{rn}.{j + 1}"
            add("causal", id=cid, ds=0, t=tidx[r[1]], g=g, orf=orf, name=(G["n"][g] or None) if g >= 0 else None,
                c=chrom, a=G["a"][g] if g >= 0 else None, e=G["e"][g] if g >= 0 else None, pcs=None, sp=sp, eff=eff,
                fdr=None, jppc=None, jv=None, q1=q, q2=-1, tier="c", pv=f"e12:S11#{rn}")
            orig[chrom]["e12:S11"][cid] = verbatim(r)
    cz = T["causal"]
    per_gene = collections.OrderedDict()
    for i in range(atlas.nlen(cz)):
        if cz["g"][i] >= 0:
            per_gene.setdefault(cz["g"][i], []).append(i)
    for g, ii in per_gene.items():
        add("cand", ds=0, g=g, tier="c", fdr=None, nt=len({cz["t"][i] for i in ii}), n=len(ii))

    # ---- dataset row ----------------------------------------------------------------------
    nq = atlas.nlen(T["qtl"])
    npk = sum(len(x) for x in T["qtl"]["xs"])
    add("datasets", id=DS_ID, sc=SC, label="Ehrenreich 2012",
        cite="Ehrenreich IM, Bloom J, Torabi N, Wang X, Jia Y, Kruglyak L (2012) Genetic Architecture of Highly Complex Chemical Resistance Traits across Four Yeast Strains. PLoS Genet 8(3): e1002570",
        doi="10.1371/journal.pgen.1002570", year=2012,
        design="X-QTL (bulk segregant analysis of millions of segregants selected for resistance) in all 6 pairwise crosses of BY, RM, YJM789 and YPS163; one selection per chemical and cross (78 experiments)",
        build="S288C reference as obtained from SGD at array design (release not stated). Positions agree with R64 gene positions at the resolution of the loci and none exceeds an R64 chromosome; used as given, no transform.",
        method="Per SNP, log10 intensity difference of MATalpha- and MATa-allele probes, minus the mean of seven cross-specific control hybridizations; Savitzky-Golay smoothing over 100-probe windows; peaks called at 1% FDR (threshold from running the caller on control data); peaks of all crosses grouped into distinct loci per trait within 200 kb of the strongest peak.",
        fdr=0.01,
        interval="None reported. Each locus is one position (the grouped locus); the grouping window was 200 kb.",
        effect="No effect sizes. Direction per cross: the parent whose allele was enriched in the resistant pool (qtl.ups, paired element by element with qtl.xs; qtl.dns holds the other parent).",
        nq=nq,
        notes=f"4 parents, 6 crosses (first-named parent MATalpha), 13 chemicals with cross-specific doses, {nq} distinct loci from {npk} per-cross peaks. "
              "Candidate genes: for the 135 allelic-singleton loci, genes with nonsynonymous SNPs in the singleton pattern (File S11), tier 'c' (candidate named by the paper); none is claimed causal.")

    atlas.save_fragment(DS_ID, dict(T))

    # ---- lazy files -------------------------------------------------------------------------
    for chrom in atlas.chrom_lengths():
        out = {"cols": {"e12:S10": qcols, "e12:S11": ccols}, "rows": {}}
        for code in ("e12:S10", "e12:S11"):
            out["rows"][code] = dict(orig[chrom].get(code, {}))
        atlas.write_json(os.path.join(atlas.DOCS, "v", f"{SC}_{chrom}.json"), out)
    extra = {}
    for code in ("e12:S1dose", "e12:S1det"):
        cols, rows = read_xlsx(code)
        extra[code] = {"cols": cols, "rows": [verbatim(r) for _, r in rows]}
    atlas.write_json(os.path.join(atlas.DOCS, "v", f"{SC}_extra.json"), extra)

    data = atlas.assemble()
    print(f"{DS_ID}: {len(pidx)} parents, {len(xidx)} crosses, {len(tidx)} traits, {nq} loci ({npk} per-cross peaks), "
          f"{atlas.nlen(T['causal'])} candidate-gene rows for {atlas.nlen(T['cand'])} genes; "
          f"ORFs not in the R64-5-1 table: {sorted(set(missing_genes))}", file=sys.stderr)
    print(f"  docs/data.json now holds {len(data['datasets']['id'])} dataset(s)", file=sys.stderr)


if __name__ == "__main__":
    main()
