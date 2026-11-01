"""Ingest Gou et al. 2019 (Genetics 211:731-740), the genetic basis of mutation rate variation in BY x RM.

    python3 -I ingest/gou2019.py          # from the repository root

Reads only from sources/: the article PDF and the unpacked supplement in sources/gou2019/raw/,
sources/gou2019/extracted/text_facts.tsv (statements quoted verbatim from the article, with page
and column, and one line of the authors' notebook), and the SGD tables in sources/sgd/.
Writes build/fragments/gou2019.json, docs/v/g19_*.json and, through ingest/atlas.py,
docs/data.json. Deterministic.

What the files hold: the four mutation-rate QTL are given in the Results text by chromosome
(XII, XIV, I, V), with the direction of effect and the variance explained; their positions,
intervals and LODs are only drawn (Figure 1B, against the genetic map) and are not tabulated
anywhere in the article, the supplement or the authors' notebooks. So the QTL rows here carry a
chromosome and no coordinates (a, e, p are null) and are listed, not drawn. The 4NQO and MMS
QTL are described as two large-effect QTL overlapping the major (chrXII, chrXIV) mutation-rate
QTL. RAD5 and MKT1 are confirmed by allele replacement (tier x) for mutation rate.

Genome build: no coordinates are given, so none applies. Gene positions shown for RAD5 and MKT1
come from the SGD R64-5-1 table.
"""
import collections
import hashlib
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import atlas  # noqa: E402

DS_ID, SC = "gou2019", "g19"
SRC = os.path.join(atlas.SOURCES, "gou2019")
RAW = os.path.join(SRC, "raw")
ART = "Gou_et_al_2019_Genetics.pdf"
SUPP = "Supplemental Figures and Tables.pdf"
S1 = "S1.File_mutation_rate_score.xlsx"
S2 = "S2.File_mutation_rate_engineer_strains.xlsx"
TF = "../extracted/text_facts.tsv"   # relative to raw/, as sources.file is


def manifest():
    out = {}
    with open(os.path.join(SRC, "MANIFEST.tsv")) as f:
        for line in f.read().splitlines()[1:]:
            name, _, sha = line.split("\t")
            out[name] = sha
    return out


def xlsx_rows(name):
    import openpyxl
    wb = openpyxl.load_workbook(os.path.join(RAW, name), read_only=True, data_only=True)
    rows = [list(r) for r in wb.worksheets[0].iter_rows(values_only=True)]
    return [("" if v is None else v) for v in rows[0]], [(i + 2, r) for i, r in enumerate(rows[1:])]


def pdf_lines(name, columns):
    """Text lines per page; for the two-column article, left then right column of each page."""
    import pdfplumber
    out = []
    with pdfplumber.open(os.path.join(RAW, name)) as pdf:
        for pi, pg in enumerate(pdf.pages):
            parts = [("L", (0, pg.width / 2)), ("R", (pg.width / 2, pg.width))] if columns else [("", (0, pg.width))]
            for side, (x0, x1) in parts:
                t = pg.crop((x0, 0, x1, pg.height)).extract_text(x_tolerance=1) or ""
                for li, l in enumerate(t.split("\n")):
                    out.append([pi + 1, side, li + 1, l])
    return out


def main():
    sha = manifest()
    sha[TF] = hashlib.sha256(open(os.path.join(SRC, "extracted", "text_facts.tsv"), "rb").read()).hexdigest()
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

    srcs = {"g19:TF": (TF, "text_facts", 1), "g19:S1": (S1, "mutation score", 1), "g19:S2": (S2, "Sheet1", 1),
            "g19:ART": (ART, "all pages, by column (text lines)", None), "g19:SUPP": (SUPP, "all pages (text lines)", None)}
    for code, (f, sheet, hrow) in srcs.items():
        add("sources", code=code, ds=0, file=f, sheet=sheet, hrow=hrow, sha256=sha[f])

    # ---- text facts -------------------------------------------------------------------------
    with open(os.path.join(SRC, "extracted", "text_facts.tsv")) as f:
        tfr = [l.rstrip("\n").split("\t") for l in f]
    tf_head = tfr[0]
    tf = [(i + 2, dict(zip(tf_head, r))) for i, r in enumerate(tfr[1:])]
    TFD = dict(tf)

    def fact(kind, contains):
        return next(n for n, d in tf if d["kind"] == kind and contains in d["quote"])

    # ---- parents, cross ---------------------------------------------------------------------
    by = add("parents", id=f"{SC}-BY", ds=0, name="BY", code="BY4724", s1="BY4724 (laboratory strain)", s2=None, p18=None,
             org="laboratory strain")
    rm = add("parents", id=f"{SC}-RM", ds=0, name="RM", code="RM11-1a", s1=None, s2="RM11-1a (vineyard strain)", p18=None,
             org="vineyard strain")
    s1_cols, s1_rows = xlsx_rows(S1)
    n_mapped = len(s1_rows)
    x = add("crosses", id=f"{SC}-BYxRM", ds=0, code="BYxRM", name="BY × RM", p1=by, p2=rm, g1=by, n=n_mapped,
            dip=None, mm=None, pv=f"g19:TF#{fact('design', '843 segregants')}")

    # ---- traits ---------------------------------------------------------------------------------
    cpb = "atlas assignment: no existing category is specific to mutation rate; 'cellular phenotype' as for Ehrenreich 2010's mitochondrial function"
    mr = add("traits", id=f"{SC}-mutation-rate", ds=0, raw="mutation rate score", name="Mutation rate (CAN1, canavanine)",
             cond=None, bt=None, ck="mutation rate", cat="cellular phenotype", catb=cpb)
    chem = "atlas assignment: DNA-damaging agent on YPD agar (Table 2); colony size as growth"
    nqo = add("traits", id=f"{SC}-4NQO", ds=0, raw="4NQO;0.1 mg/ml", name="4-NQO 0.1 mg/mL", cond="0.1 mg/ml", bt=None,
              ck="4nqo 0.1mg/mL", cat="drug/chemical", catb=chem)
    mms = add("traits", id=f"{SC}-MMS", ds=0, raw="MMS;0.01%", name="MMS 0.01%", cond="0.01%", bt=None,
              ck="mms 0.01%", cat="drug/chemical", catb=chem)
    add("traits", id=f"{SC}-H2O2", ds=0, raw="H2O2;4 mM", name="Hydrogen peroxide 4 mM", cond="4 mM", bt=None,
        ck="hydrogen peroxide 4mM", cat="drug/chemical", catb=chem)

    G = atlas.genes()
    gname = {}
    for gi, n in enumerate(G["n"]):
        if n:
            gname.setdefault(n, gi)

    # ---- QTL ----------------------------------------------------------------------------------
    # Mutation rate: chromosome, direction and variance explained from the Results text.
    ve_row = fact("ve", "8.8, 6.1, 3.1")
    rm_row, by_row = fact("direction", "XII and I, the RM allele"), fact("direction", "XIV and V, the BY allele")
    MR = [("chrXII", 0.088, rm, "RAD5"), ("chrXIV", 0.061, by, "MKT1"), ("chrI", 0.031, rm, None), ("chrV", 0.026, by, None)]
    orig = collections.defaultdict(lambda: collections.defaultdict(dict))
    causal = []
    for k, (chrom, ve, up, gene) in enumerate(MR):
        qid = f"{SC}-w-{k + 1}"
        g = gname[gene] if gene else -1
        q = add("qtl", id=qid, ds=0, t=mr, x=x, m="w", c=chrom, a=None, e=None, p=None, cm=None, lod=None, q=None, ef=None,
                ve=ve, lv=None, af=None, cc=None, pa=None, anc=None, daf=None, up=up, ups=None, dns=None, xs=None,
                k="mutation-rate QTL (position drawn in Figure 1B only)", g=g, gf=None, gt="x" if gene else None,
                pv=f"g19:TF#{ve_row}")
        orig[chrom]["g19:TF"][qid] = [[TFD[r][h] for h in tf_head] for r in (ve_row, rm_row if up == rm else by_row)]
        if gene:
            causal.append((mr, g, q, fact("candidate", f"{'RM RAD5' if gene == 'RAD5' else 'BY MKT1'} allele"), "x"))
    # 4NQO and MMS: two large-effect QTL overlapping the major (chrXII, chrXIV) mutation-rate QTL.
    mut_row = fact("mutagen", "two large-effect QTL for 4NQO and MMS")
    nqo_genes = fact("mutagen", "causal genes for 4NQO sensitivity")
    n = len(MR)
    for t, tname in ((nqo, "4NQO"), (mms, "MMS")):
        for chrom, gene in (("chrXII", "RAD5"), ("chrXIV", "MKT1")):
            n += 1
            qid = f"{SC}-w-{n}"
            g = gname[gene] if t == nqo else -1
            q = add("qtl", id=qid, ds=0, t=t, x=x, m="w", c=chrom, a=None, e=None, p=None, cm=None, lod=None, q=None, ef=None,
                    ve=None, lv=None, af=None, cc=None, pa=None, anc=None, daf=None, up=-1, ups=None, dns=None, xs=None,
                    k=f"{tname} QTL overlapping the {chrom} mutation-rate QTL (text; drawn in Figure 4 only)",
                    g=g, gf=None, gt="c" if t == nqo else None, pv=f"g19:TF#{mut_row}")
            orig[chrom]["g19:TF"][qid] = [[TFD[mut_row][h] for h in tf_head]]
            if t == nqo:
                causal.append((nqo, g, q, nqo_genes, "c"))

    for t, g, q, row, tier in causal:
        cid = f"{SC}-c-{t}-{G['id'][g]}"
        add("causal", id=cid, ds=0, t=t, g=g, orf=G["id"][g], name=G["n"][g], c=G["c"][g], a=G["a"][g], e=G["e"][g],
            pcs=None, sp=None, eff=None, fdr=None, jppc=None, jv=None, q1=q, q2=-1, tier=tier, pv=f"g19:TF#{row}")
        orig[G["c"][g]]["g19:TF"][cid] = [[TFD[row][h] for h in tf_head]]
    cz = T["causal"]
    rank = {"x": 0, "v": 1, "5": 2, "c": 3, "20": 4}
    per_gene = C()
    for i in range(atlas.nlen(cz)):
        per_gene.setdefault(cz["g"][i], []).append(i)
    for g, ii in per_gene.items():
        add("cand", ds=0, g=g, tier=min((cz["tier"][i] for i in ii), key=rank.get), fdr=None,
            nt=len({cz["t"][i] for i in ii}), n=len(ii))

    # ---- variance explained (text) -----------------------------------------------------------
    add("vexp", t=mr, x=x, wve=0.207, jve=None, h2=0.30, h2se=None, re=None, rese=None,
        pv=f"g19:TF#{fact('ve_total', '20.7%')}")

    add("datasets", id=DS_ID, sc=SC, label="Gou 2019",
        cite="Gou L, Bloom JS, Kruglyak L. 2019. The Genetic Basis of Mutation Rate Variation in Yeast. Genetics 211:731-740",
        doi="10.1534/genetics.118.301609", year=2019,
        design="segregant panel, BY × RM (1040 segregants genotyped in Bloom et al. 2013; 843 mapped); spontaneous mutation rate "
               "at CAN1 by a high-throughput fluctuation assay, and growth on 4NQO, MMS and H2O2",
        build="none applies: the files give no QTL coordinates (chromosome only). Gene positions for RAD5 and MKT1 are from SGD R64-5-1.",
        method="Single-marker LOD from the Pearson correlation between genotype and the mutation rate score (residual of m after "
               "cells per culture and plate) or growth residual; 5% FWER thresholds from 1000 permutations: 3.52 (mutation rate), "
               "3.62 (4NQO), 3.61 (MMS), 3.64 (H2O2). QTL-QTL interactions by an F test of an interaction term.",
        fdr=None,
        interval="1.5-LOD drop in the paper, but no interval is tabulated: positions are drawn only (Figure 1B, Figure 4). "
                 "QTL rows here have a chromosome and no coordinates.",
        effect="Not tabulated. Direction from the text: on chrXII and chrI the RM allele, on chrXIV and chrV the BY allele confers a "
               "higher mutation rate (qtl.up). Variance explained per QTL from the text (qtl.ve).",
        nq=atlas.nlen(T["qtl"]),
        notes=f"1 cross (BY × RM), {n_mapped} segregants mapped of 1040, 4 traits (mutation rate; 4NQO, MMS, H2O2 growth). "
              "Four mutation-rate QTL (chrXII 8.8%, chrXIV 6.1%, chrI 3.1%, chrV 2.6%; 20.7% together, 69% of the 30% additive "
              "heritability); RAD5 and MKT1 confirmed by allele replacement. 4NQO and MMS QTL overlap chrXII and chrXIV. "
              "Positions are only plotted, so these QTL are listed but not drawn on the genome.")

    atlas.save_fragment(DS_ID, T)

    # ---- lazy files ------------------------------------------------------------------------------
    for chrom in atlas.chrom_lengths():
        rows = {}
        for rid, vals in orig[chrom].get("g19:TF", {}).items():
            rows[rid] = vals[0]   # the row the provenance key points to
        atlas.write_json(os.path.join(atlas.DOCS, "v", f"{SC}_{chrom}.json"),
                         {"cols": {"g19:TF": tf_head}, "rows": {"g19:TF": rows}})
    s2_cols, s2_rows = xlsx_rows(S2)
    extra = {"g19:TF": {"cols": tf_head, "rows": [[d[h] for h in tf_head] for _, d in tf]},
             "g19:S1": {"cols": s1_cols, "rows": [r for _, r in s1_rows]},
             "g19:S2": {"cols": s2_cols, "rows": [r for _, r in s2_rows]},
             "g19:ART": {"cols": ["page", "column", "line", "text"], "rows": pdf_lines(ART, True)},
             "g19:SUPP": {"cols": ["page", "column", "line", "text"], "rows": pdf_lines(SUPP, False)}}
    atlas.write_json(os.path.join(atlas.DOCS, "v", f"{SC}_extra.json"), extra)

    data = atlas.assemble()
    print(f"{DS_ID}: 1 cross ({n_mapped} segregants mapped), {atlas.nlen(T['traits'])} traits, {atlas.nlen(T['qtl'])} QTL "
          f"(chromosome only), {atlas.nlen(T['causal'])} candidate rows ({atlas.nlen(T['cand'])} genes); "
          f"data.json has {len(data['datasets']['id'])} dataset(s)", file=sys.stderr)


if __name__ == "__main__":
    main()
