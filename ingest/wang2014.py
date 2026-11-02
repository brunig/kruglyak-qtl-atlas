"""Ingest Wang X & Kruglyak L 2014 (PLoS Genet 10:e1004894), haloperidol resistance in BY x RM.

    python3 -I ingest/wang2014.py          # from the repository root

Reads only from sources/: the supporting-information files in sources/wang2014/raw/, their
verbatim flattenings in sources/wang2014/extracted/S*_Table.tsv (made by extract_tables.py and
checked against the .docx files by the tests), the quoted sentences in extracted/text_facts.tsv,
and the SGD tables in sources/sgd/. Writes build/fragments/wang2014.json, docs/v/w14_*.json and,
through ingest/atlas.py, docs/data.json. Deterministic.

Design: linkage in 1008 individually genotyped BY x RM segregants (the Bloom et al. 2013 panel),
growth on five haloperidol doses; R/qtl (Haley-Knott and non-parametric), QTL called at LOD 3
with 1000 permutations; nine loci. So m = "w" (not X-QTL). S1 Table gives, per dose, every locus
in the additive model with its drop-one LOD and % variance explained; loci are named by
chromosome and arm (1L, 1R, 2, 7, 12, 13, 14, 15-1, 15-2). No position, peak marker or interval
is tabulated anywhere (Figure 3A draws the 1.5-LOD intervals), so QTL rows have a chromosome and
no coordinates (a, e, p null): listed, not drawn.

Direction (text): the RM allele raises resistance at the chrI right-arm locus and at "loci on
chromosomes V, XII and XV"; BY alleles at "the remaining five loci". S1 Table has no chrV locus.
Most defensible reading (OPEN_QUESTIONS): RM at 1R, 12, 15-1 and 15-2 (four loci), BY at 1L, 2,
7, 13 and 14 (five, matching "the remaining five"). SWH1 (chrI right arm), MKT1 (chrXIV) and IRA2
(chrXV) are confirmed by allele replacement (tier x); IRA2 is assigned to locus 15-1, the chrXV
locus that is major at the higher doses as the text says of the chrXV QTL (reading, labelled).

Genome build: no QTL coordinates. S2 Table's SWH1 SNP positions fall inside YAR042W on R64-5-1,
so the paper's coordinate system for that table is R64-compatible; nothing is transformed.
"""
import collections
import csv
import hashlib
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import atlas  # noqa: E402

DS_ID, SC = "wang2014", "w14"
SRC = os.path.join(atlas.SOURCES, "wang2014")
RAW = os.path.join(SRC, "raw")
EXT = os.path.join(SRC, "extracted")
TABLES = {  # code -> (docx in raw/, flattened TSV in extracted/, caption)
    "w14:S1": ("journal.pgen.1004894.s002.docx", "S1_Table.tsv", "S1 Table: additive QTL model per dose"),
    "w14:S2": ("journal.pgen.1004894.s003.docx", "S2_Table.tsv", "S2 Table: SWH1 coding polymorphisms"),
    "w14:S3": ("journal.pgen.1004894.s004.docx", "S3_Table.tsv", "S3 Table: QTL models with two-way interactions"),
    "w14:S4": ("journal.pgen.1004894.s005.docx", "S4_Table.tsv", "S4 Table: allele-replacement ANOVA, both backgrounds"),
    "w14:S5": ("journal.pgen.1004894.s006.docx", "S5_Table.tsv", "S5 Table: allele-replacement ANOVA, BY background"),
    "w14:S6": ("journal.pgen.1004894.s007.docx", "S6_Table.tsv", "S6 Table: Tukey-HSD tests"),
}
TF = "../extracted/text_facts.tsv"   # relative to raw/, as sources.file is
DOSES = ["40", "80", "120", "160", "200"]
LOCUS = {  # S1/S3 label -> (chromosome, description)
    "1L": ("chrI", "chrI left arm"), "1R": ("chrI", "chrI right arm"), "2": ("chrII", "chrII"),
    "7": ("chrVII", "chrVII"), "12": ("chrXII", "chrXII"), "13": ("chrXIII", "chrXIII"), "14": ("chrXIV", "chrXIV"),
    "15-1": ("chrXV", "chrXV locus 15-1"), "15-2": ("chrXV", "chrXV locus 15-2"),
}
RM_UP = {"1R", "12", "15-1", "15-2"}       # reading of the direction sentence; see module docstring
GENE_OF = {"1R": "SWH1", "14": "MKT1", "15-1": "IRA2"}


def manifest():
    out = {}
    with open(os.path.join(SRC, "MANIFEST.tsv")) as f:
        for line in f.read().splitlines()[1:]:
            name, _, sha = line.split("\t")
            out[name] = sha
    return out


def tsv(path):
    with open(path, encoding="utf-8") as f:
        rows = list(csv.reader(f, delimiter="\t", quoting=csv.QUOTE_NONE))
    return rows[0], [(i + 2, r) for i, r in enumerate(rows[1:])]   # (line number in the TSV, header = 1)


def dose_of(heading):
    return heading.split()[0]   # "40 μM haloperidol:" -> "40"


def fnum(v):
    v = (v or "").strip()
    if not v:
        return None
    try:
        return float(v)
    except ValueError:
        return None


def main():
    sha = manifest()
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

    for code, (docx, flat, cap) in TABLES.items():
        add("sources", code=code, ds=0, file=docx, sheet=f"{cap} (Word tables, flattened verbatim in extracted/{flat})", hrow=1,
            sha256=sha[docx])
    add("sources", code="w14:TF", ds=0, file=TF, sheet="sentences quoted from article.html", hrow=1,
        sha256=hashlib.sha256(open(os.path.join(EXT, "text_facts.tsv"), "rb").read()).hexdigest())

    tf_head, tf = tsv(os.path.join(EXT, "text_facts.tsv"))

    def fact(kind, contains):
        return next(n for n, r in tf if r[0] == kind and contains in r[1])

    # ---- parents, cross ---------------------------------------------------------------------
    design = fact("design", "1008 segregants from a cross")
    by = add("parents", id=f"{SC}-BY", ds=0, name="BY", code="BY4724", s1="BY4724 (laboratory strain)", s2=None, p18=None,
             org="laboratory strain")
    rm = add("parents", id=f"{SC}-RM", ds=0, name="RM", code="RM11-1a", s1=None, s2="RM11-1a (vineyard strain)", p18=None,
             org="vineyard strain")
    x = add("crosses", id=f"{SC}-BYxRM", ds=0, code="BYxRM", name="BY × RM", p1=by, p2=rm, g1=None, n=1008,
            dip=None, mm=None, pv=f"w14:TF#{design}")

    # ---- traits: one per dose -----------------------------------------------------------------
    catb = "atlas assignment: drug added to YPD agar (Methods); growth relative to a DMSO control"
    tidx = {}
    for d in DOSES:
        tidx[d] = add("traits", id=f"{SC}-haloperidol-{d}uM", ds=0, raw=f"haloperidol;{d} µM", name=f"Haloperidol {d} µM",
                      cond=f"{d} µM", bt=None, ck=f"haloperidol {d}uM", cat="drug/chemical", catb=catb)

    G = atlas.genes()
    gname = {}
    for gi, n in enumerate(G["n"]):
        if n:
            gname.setdefault(n, gi)

    # ---- QTL: S1 Table, one row per (dose, locus) -----------------------------------------------
    s1_head, s1 = tsv(os.path.join(EXT, "S1_Table.tsv"))
    orig = collections.defaultdict(lambda: collections.defaultdict(dict))
    qkey = {}
    dir_rm, dir_1r = fact("direction", "chromosomes V, XII and XV"), fact("direction", "right arm of chromosome I")
    for line, r in s1:
        label = r[4]
        if label in ("", "Total") or r[5] == "df":
            continue
        d = dose_of(r[3])
        chrom, desc = LOCUS[label]
        gene = GENE_OF.get(label)
        qid = f"{SC}-w-{line}"
        q = add("qtl", id=qid, ds=0, t=tidx[d], x=x, m="w", c=chrom, a=None, e=None, p=None, cm=None,
                lod=atlas.r6(fnum(r[7])), q=None, ef=None, ve=atlas.r6(fnum(r[8]) / 100), lv=None, af=None, cc=None, pa=None,
                anc=None, daf=None, up=rm if label in RM_UP else by, ups=None, dns=None, xs=None,
                k=f"locus {label} ({desc}); position drawn in Figure 3A only",
                g=gname[gene] if gene else -1, gf=None, gt="x" if gene else None, pv=f"w14:S1#{line}")
        qkey[(d, label)] = q
        orig[chrom]["w14:S1"][qid] = r

    # ---- candidate genes: confirmed by allele replacement ----------------------------------------
    gene_fact = {"SWH1": fact("candidate", "SWH1 RM is the resistant allele"),
                 "MKT1": fact("candidate", "we confirmed MKT1"),
                 "IRA2": fact("candidate", "The IRA2 RM allele")}
    for (d, label), q in qkey.items():
        gene = GENE_OF.get(label)
        if not gene:
            continue
        g = gname[gene]
        cid = f"{SC}-c-{d}-{gene}"
        add("causal", id=cid, ds=0, t=tidx[d], g=g, orf=G["id"][g], name=gene, c=G["c"][g], a=G["a"][g], e=G["e"][g],
            pcs=None, sp=None, eff=None, fdr=None, jppc=None, jv=None, q1=q, q2=-1, tier="x", pv=f"w14:TF#{gene_fact[gene]}")
        orig[G["c"][g]]["w14:TF"][cid] = tf[gene_fact[gene] - 2][1]
    cz = T["causal"]
    per_gene = C()
    for i in range(atlas.nlen(cz)):
        per_gene.setdefault(cz["g"][i], []).append(i)
    for g, ii in per_gene.items():
        add("cand", ds=0, g=g, tier="x", fdr=None, nt=len({cz["t"][i] for i in ii}), n=len(ii))

    # ---- variance explained per dose: the S1 "Total" rows ---------------------------------------
    for line, r in s1:
        if r[4] == "Total":
            add("vexp", t=tidx[dose_of(r[3])], x=x, wve=atlas.r6(fnum(r[8]) / 100), jve=None, h2=None, h2se=None, re=None,
                rese=None, pv=f"w14:S1#{line}")

    nq = atlas.nlen(T["qtl"])
    add("datasets", id=DS_ID, sc=SC, label="Wang 2014",
        cite="Wang X, Kruglyak L. 2014. Genetic basis of haloperidol resistance in Saccharomyces cerevisiae is complex and dose dependent. PLoS Genet 10:e1004894",
        doi="10.1371/journal.pgen.1004894", year=2014,
        design="segregant panel, BY4724 × RM11-1a, 1008 segregants (the Bloom et al. 2013 panel); end-point colony growth on YPD "
               "with 40, 80, 120, 160 or 200 µM haloperidol, normalized to a DMSO control",
        build="none applies to the QTL: no coordinates are tabulated (chromosome only). S2 Table's SWH1 SNP positions lie inside "
              "YAR042W on R64-5-1, so that table is R64-compatible; nothing is transformed.",
        method="Linkage by R/qtl (Haley-Knott regression and non-parametric mapping) on control-normalized residuals; QTL called at "
               "LOD 3, significance by 1000 permutations; per-dose additive model (S1 Table: drop-one LOD and % variance); "
               "two-way interactions at Bonferroni p < 0.005 (S3 Table).",
        fdr=None,
        interval="1.5-LOD drop in the paper (Figure 3A), but no interval or position is tabulated: QTL rows carry a chromosome "
                 "(and arm, in qtl.k) and no coordinates.",
        effect="Not tabulated (Figure 3B plots effects; positive = the BY allele raises growth in haloperidol). Direction per "
               "locus from the text in qtl.up: RM at chrI right arm, chrXII and both chrXV loci; BY at the other five "
               "(reading; the text also names a chrV locus that S1 Table does not have). LOD and variance are the drop-one "
               "values of the per-dose additive model.",
        nq=nq,
        notes="1 cross (BY × RM), 1008 segregants, 5 haloperidol doses, 9 loci (chrI left and right arms, II, VII, XII, XIII, XIV, "
              "and two on XV), 32 locus-dose rows from S1 Table. SWH1 (chrI right arm), MKT1 (chrXIV) and IRA2 (chrXV) confirmed "
              "by allele replacement. Positions are only plotted, so these QTL are listed but not drawn on the genome.")

    atlas.save_fragment(DS_ID, T)

    # ---- lazy files -------------------------------------------------------------------------------
    for chrom in atlas.chrom_lengths():
        cols, rows = {}, {}
        for code in ("w14:S1", "w14:TF"):
            cols[code] = s1_head if code == "w14:S1" else tf_head
            rows[code] = dict(orig[chrom].get(code, {}))
        atlas.write_json(os.path.join(atlas.DOCS, "v", f"{SC}_{chrom}.json"), {"cols": cols, "rows": rows})
    extra = {"w14:TF": {"cols": tf_head, "rows": [r for _, r in tf]}}
    for code, (_, flat, _) in TABLES.items():
        h, rs = tsv(os.path.join(EXT, flat))
        extra[code] = {"cols": h, "rows": [r for _, r in rs]}
    atlas.write_json(os.path.join(atlas.DOCS, "v", f"{SC}_extra.json"), extra)

    data = atlas.assemble()
    print(f"{DS_ID}: 1 cross (1008 segregants), {len(DOSES)} traits, {nq} QTL rows (chromosome only, "
          f"{len({k[1] for k in qkey})} loci), {atlas.nlen(T['causal'])} candidate rows; "
          f"data.json has {len(data['datasets']['id'])} dataset(s)", file=sys.stderr)


if __name__ == "__main__":
    main()
